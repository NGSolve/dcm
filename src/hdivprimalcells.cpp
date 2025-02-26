#include <comp.hpp>
#include "hdivprimalcells.hpp"
#include "intrules.hpp"
#include "dcs_common.hpp"
#include "supersparse.hpp"
#include "smallsparse.hpp"

namespace ngcomp
{

 class ApplyMassHDivPrimalCells : public ApplyMass
  {
    Matrix<> eval; 
    SmallSparseMatrix speval, spevalT;
    bool use_sparseeval = false;
    IntegrationRule trigir;
    Matrix<Mat<2,2>> rho_eval;
    public:
      ApplyMassHDivPrimalCells (shared_ptr<FESpace> afes,
          shared_ptr<CoefficientFunction> arho,
          bool ainverse,
          shared_ptr<Region> adefinedon,
          LocalHeap & alh);

      ~ApplyMassHDivPrimalCells();

      void Mult(const BaseVector & vec, BaseVector &prod) const;

      void MultAdd(double val, const BaseVector & vec, BaseVector & prod) const;

      shared_ptr<BaseMatrix> InverseMatrix(shared_ptr<BitArray> subset) const 
      {
        return make_shared<ApplyMassHDivPrimalCells> (fes, rho, !inverse, definedon, lh);
      }
  };

  
  ApplyMassHDivPrimalCells::
      ApplyMassHDivPrimalCells (shared_ptr<FESpace> afes,
          shared_ptr<CoefficientFunction> arho,
          bool ainverse,
          shared_ptr<Region> adefinedon,
          LocalHeap & alh)
        : ApplyMass(afes, arho, ainverse, adefinedon, alh)
      { 
        cout << IM(3) << "ApplyMassHDivPrimalCells::ApplyMassHDivPrimalCells called" << endl; 
        cout << IM(3) << " inverse = " << inverse << endl;

      auto feshdiv = dynamic_pointer_cast<HDivPrimalCells>(fes); 
      auto irs = feshdiv->GetIntegrationRules();
      trigir = std::move(irs[ET_TRIG]);
      //cout << trigir << endl;

      auto & felref = static_cast<const HDivPrimalCellTrig&>(fes->GetFE(ElementId(0), lh));
      Matrix<> shaperef(felref.GetNDof(), 2);
      
      Matrix<> shapes(felref.GetNDof(), 2*trigir.Size());
      for (int i = 0; i < trigir.Size(); i++)
        felref.CalcShape (trigir[i], shapes.Cols(2*i, 2*i+2));


      if (inverse)
      { 
        Matrix<> matref(felref.GetNDof());
        matref = 0.0;
        for (int i = 0; i < trigir.Size(); i++)
          {
            auto shaperef = shapes.Cols(2*i, 2*i+2);
            matref += trigir[i].Weight() * shaperef * Trans(shaperef);
          }
        CalcInverse (matref);
        
        eval = Matrix(Trans(shapes)*matref);
      }
      else
        eval = Matrix(Trans(shapes));

      for (int i : Range(eval.Height()))
        for (int j : Range(eval.Width()))
          if (fabs(eval(i,j)) < 1e-10)
            eval(i,j) = 0;

      Array<Mat<2,2>> piola_factors(trigir.Size());
      if (inverse)
        {
          // cout << "eval = " << endl << eval << endl;
          for (int i = 0; i < trigir.Size(); i++)
            {
              auto piola = felref.GetPiola (trigir[i]);
              Mat<2,2> fac = trigir[i].Weight()*Trans(piola);
              eval.Rows(2*i, 2*i+2) = Matrix(fac * eval.Rows(2*i, 2*i+2 ));
              piola_factors[i] = fac;
            }

          for (int i : Range(eval.Height()))
            for (int j : Range(eval.Width()))
              if (fabs(eval(i,j)) < 1e-10)
                eval(i,j) = 0;
          cout << IM(3) << "eval, now = " << endl << eval << endl;          
        }


      
      speval = SmallSparseMatrix(eval, 1e-8);
      spevalT = SmallSparseMatrix(Matrix(Trans(eval)), 1e-8);
      cout << IM(3) << "sparse-eval nze = " << speval.NZE() << " / " << eval.Height() * eval.Width() << endl;
      use_sparseeval = 10*speval.NZE() < eval.Height() * eval.Width();
      
      // cout << "speval = " << speval << endl;

      rho_eval = Matrix(fes->GetMeshAccess()->GetNE(),trigir.Size());
      IterateElements(*fes, VOL, alh,
          [&](FESpace::Element el, LocalHeap & lh)
          {
             //HeapReset hr(lh);
             const ElementTransformation & trafo = el.GetTrafo();
             MappedIntegrationRule<2,2> mir(trigir, trafo, lh);
             auto rhois = rho_eval.Row(el.Nr());
             for (int i = 0; i < trigir.Size(); i++)
               {
                 Mat<2,2> rhoi = Identity(2);
                 if (rho)
                 {
                  if (rho->Dimension() == 1)
                    rhoi*=rho->Evaluate(mir[i]);
                  else if (rho->Dimension() == 4) //or trans(rho)???
                    rho->Evaluate(mir[i],FlatVector<> (4,&rhoi(0,0)));
                  else 
                  {
                  ostringstream ost;
                  ost << "wrong dimension of coefficient (dim = "<< rho->Dimension()<< "), allowed are 1 or 4"<< endl;
                  throw Exception (ost.str());
                  }
                 }
                 Mat<2,2> trans = mir[i].GetJacobian();
                 
                 rhoi = 1/mir[i].GetJacobiDet()*Trans(trans)*rhoi*trans;
                 if (inverse)
                   {
                     rhoi = piola_factors[i] * rhoi * Trans(piola_factors[i]);
                     rhoi = Inv(rhoi);
                   }
                 rhoi *= trigir[i].Weight();
                 //cout << "i = " << i << "mir[i].GetWeight() = "<< mir[i].GetWeight() << "ir[i].w*det(trans)" << trigir[i].Weight()*fabs(Det(trans))<< endl;
                rhois(i)=rhoi;
               }
          
          });
      //cout << "rho_vals: " << rho_eval << endl;
      
    }
    ApplyMassHDivPrimalCells::~ApplyMassHDivPrimalCells() { ; }

    void ApplyMassHDivPrimalCells::Mult(const BaseVector & vec, BaseVector &prod) const
    {
      prod = 0;
      MultAdd(1.,vec,prod);
      }
    void ApplyMassHDivPrimalCells::MultAdd(double val, const BaseVector & vec, BaseVector &prod) const
    { 
      static Timer t("ApplyMassHDiv::MultAdd"); RegionTimer reg(t);
      static Timer tmult("ApplyMassHDiv::MultAdd mult");
      static Timer tmultT("ApplyMassHDiv::MultAdd multtrans");
      //cout << "Multadd called with val = " << val << endl;
      auto ma = fes->GetMeshAccess(); 
      auto vecmat = vec.FV<double>().AsMatrix(ma->GetNE(), eval.Width()); 
      auto prodmat = prod.FV<double>().AsMatrix(ma->GetNE(), eval.Width()); 

      ParallelForRange
        (ma->GetNE(), [&](IntRange myrange)
         {
           constexpr size_t BS = 96;
           Matrix tmpmat(BS, eval.Height());

           for (size_t first = myrange.First(); first < myrange.Next(); first += BS)
             {
               size_t next = min(first+BS, myrange.Next());
               auto r = Range(first, next);
               auto allpntvals = tmpmat.Rows(r.Size());
               
               {
                 RegionTracer tr(TaskManager::GetThreadId(), tmult);             
                 if (use_sparseeval)
                   speval.Mult (Trans(vecmat.Rows(r)), Trans(allpntvals));
                 else
                   allpntvals = vecmat.Rows(r) * Trans(eval);
               }
               
               for (auto elnr : r)
                 {
                   auto pntvals = allpntvals.Row(elnr-r.First());
                   auto rhois = rho_eval.Row(elnr);
                   
                   for (int i = 0; i < trigir.Size(); i++)
                     {
                       Vec<2> rhopval = rhois(i) * pntvals.Range(2*i,2*i+2);
                       pntvals.Range(2*i,2*i+2) = rhopval;
                     }
                 }

               {
                 RegionTracer tr(TaskManager::GetThreadId(), tmultT);
                 if (use_sparseeval)
                   spevalT.MultAdd (val, Trans(allpntvals), Trans(prodmat.Rows(r)));
                 else
                   {
                     allpntvals *= val;
                     prodmat.Rows(r) += allpntvals * eval;
                   }
               }
             }
         },
         TasksPerThread(4));
    }
  
  Mat<2,2> HDivPrimalCellTrig:: GetPiola (const IntegrationPoint & ip) const
  {
    int ndof_edge = IR.Size();
    int ndof_quad = 2*IR.Size()*(IR.Size()-1);
    
    double lam[] = { ip(0), ip(1), 1-ip(0)-ip(1) };
    int maxlam = 0;
    for (int j = 1; j < 3; j++)
      if (lam[j] > lam[maxlam])
        maxlam = j;
    
    double lx = lam[(maxlam+1)%3];
    double ly = lam[(maxlam+2)%3];
    double xi  = (3+2*lx-2*ly)/2 - sqrt( sqr( (2*ly-2*lx-3)/2 ) - 6*lx );      
    double eta = (3+2*ly-2*lx)/2 - sqrt( sqr( (2*lx-2*ly-3)/2 ) - 6*ly );
    
    
    Vec<2> verts[3] = { { 1, 0 }, { 0, 1 }, { 0, 0 } };      
    Vec<2> p0 = verts[maxlam];
    Vec<2> p1 = 0.5 * (verts[maxlam]+verts[(maxlam+1)%3]);
    Vec<2> p2 = { 1.0/3, 1.0/3 };
    Vec<2> p3 = 0.5 * (verts[maxlam]+verts[(maxlam+2)%3]);
    
    Vec<2> dpdxi = eta*(p2-p3)+(1-eta)*(p1-p0);
    Vec<2> dpdeta = xi*(p2-p1)+(1-xi)*(p3-p0);
    
    Mat<2,2> dxy_dxieta;
    dxy_dxieta.Col(0) = dpdxi;
    dxy_dxieta.Col(1) = dpdeta;
    Mat<2,2> trafo = 1.0/Det(dxy_dxieta)*dxy_dxieta; // Piola
    
    return trafo;
  }
  
  double HDivPrimalCellTrig:: GetDet (const IntegrationPoint & ip) const
  {
    int ndof_edge = IR.Size();
    int ndof_quad = 2*IR.Size()*(IR.Size()-1);
    
    double lam[] = { ip(0), ip(1), 1-ip(0)-ip(1) };
    int maxlam = 0;
    for (int j = 1; j < 3; j++)
      if (lam[j] > lam[maxlam])
        maxlam = j;
    
    double lx = lam[(maxlam+1)%3];
    double ly = lam[(maxlam+2)%3];
    double xi  = (3+2*lx-2*ly)/2 - sqrt( sqr( (2*ly-2*lx-3)/2 ) - 6*lx );      
    double eta = (3+2*ly-2*lx)/2 - sqrt( sqr( (2*lx-2*ly-3)/2 ) - 6*ly );
    
    
    Vec<2> verts[3] = { { 1, 0 }, { 0, 1 }, { 0, 0 } };      
    Vec<2> p0 = verts[maxlam];
    Vec<2> p1 = 0.5 * (verts[maxlam]+verts[(maxlam+1)%3]);
    Vec<2> p2 = { 1.0/3, 1.0/3 };
    Vec<2> p3 = 0.5 * (verts[maxlam]+verts[(maxlam+2)%3]);
    
    Vec<2> dpdxi = eta*(p2-p3)+(1-eta)*(p1-p0);
    Vec<2> dpdeta = xi*(p2-p1)+(1-xi)*(p3-p0);
    
    Mat<2,2> dxy_dxieta;
    dxy_dxieta.Col(0) = dpdxi;
    dxy_dxieta.Col(1) = dpdeta;
    return Det(dxy_dxieta);
  }
  
    void HDivPrimalCellTrig:: CalcShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<> shape) const
    {
      // int ndof_edge = order+1;
      // int ndof_quad = 2*sqr(order+1);
      int ndof_edge = IR.Size();
      int ndof_quad = 2*IR.Size()*(IR.Size()-1);
      
      ArrayMem<double, 20> polx(IR.Size()), poly(IR.Size());      
      
      double lam[] = { ip(0), ip(1), 1-ip(0)-ip(1) };
      int maxlam = 0;
      for (int j = 1; j < 3; j++)
        if (lam[j] > lam[maxlam])
          maxlam = j;
      
      //shape = 0;
      shape.AddSize(ndof, 2) = 0;

      double lx = lam[(maxlam+1)%3];
      double ly = lam[(maxlam+2)%3];
      double xi  = (3+2*lx-2*ly)/2 - sqrt( sqr( (2*ly-2*lx-3)/2 ) - 6*lx );      
      double eta = (3+2*ly-2*lx)/2 - sqrt( sqr( (2*lx-2*ly-3)/2 ) - 6*ly );

      LagrangePolynomials(xi, IR, polx);
      LagrangePolynomials(eta, IR, poly);

      Vec<2> verts[3] = { { 1, 0 }, { 0, 1 }, { 0, 0 } };      
      Vec<2> p0 = verts[maxlam];
      Vec<2> p1 = 0.5 * (verts[maxlam]+verts[(maxlam+1)%3]);
      Vec<2> p2 = { 1.0/3, 1.0/3 };
      Vec<2> p3 = 0.5 * (verts[maxlam]+verts[(maxlam+2)%3]);

      Vec<2> dpdxi = eta*(p2-p3)+(1-eta)*(p1-p0);
      Vec<2> dpdeta = xi*(p2-p1)+(1-xi)*(p3-p0);

      Mat<2,2> dxy_dxieta;
      dxy_dxieta.Col(0) = dpdxi;
      dxy_dxieta.Col(1) = dpdeta;
      Mat<2,2> trafo = 1.0/Det(dxy_dxieta)*dxy_dxieta; // Piola
      
      auto assign_shapex = [trafo, shape, &polx, &poly](int nr, int ix, int iy)
        {
          shape.Row(nr) = trafo*Vec<2>(polx[ix]*poly[iy], 0);
        };
      auto assign_shapey = [trafo, shape, &polx, &poly](int nr, int ix, int iy)
        {
          shape.Row(nr) = trafo*Vec<2>(0, -polx[ix]*poly[iy]);
        };

      int ii = 0;

      // x-edge
      ii = (maxlam+2)%3 * ndof_edge;
      for (int i = 0; i < poly.Size(); i++)      
        assign_shapex(ii++, 0, i);
      
      ii = (maxlam+1)%3 * ndof_edge;
      for (int i = 0; i < polx.Size(); i++)
        assign_shapey(ii++, i, 0);
      
      ii = 3*ndof_edge + maxlam*ndof_quad;
      for (int i = 1; i < polx.Size(); i++)
        for (int j = 0; j < poly.Size(); j++)
          assign_shapex(ii++, i, j);
      for (int j = 1; j < poly.Size(); j++)
        for (int i = 0; i < polx.Size(); i++)
          assign_shapey(ii++, i, j);
    }
    void HDivPrimalCellTrig:: CalcAltShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<> shape) const
    {

      //
      // int ndof_edge = order+1;
      // int ndof_quad = 2*sqr(order+1);
      int ndof_edge = IR.Size();
      int ndof_quad = 2*IR.Size()*(IR.Size()-1);
      
      ArrayMem<double, 20> polx(IR.Size()), poly(IR.Size());      
      
      double lam[] = { ip(0), ip(1), 1-ip(0)-ip(1) };
      int maxlam = 0;
      for (int j = 1; j < 3; j++)
        if (lam[j] > lam[maxlam])
          maxlam = j;
      
      //shape = 0;
      shape.AddSize(ndof, 2) = 0;

      double lx = lam[(maxlam+1)%3];
      double ly = lam[(maxlam+2)%3];
      double xi  = (3+2*lx-2*ly)/2 - sqrt( sqr( (2*ly-2*lx-3)/2 ) - 6*lx );      
      double eta = (3+2*ly-2*lx)/2 - sqrt( sqr( (2*lx-2*ly-3)/2 ) - 6*ly );

      LagrangePolynomials(xi, IR, polx);
      LagrangePolynomials(eta, IR, poly);

      Vec<2> verts[3] = { { 1, 0 }, { 0, 1 }, { 0, 0 } };      
      Vec<2> p0 = verts[maxlam];
      Vec<2> p1 = 0.5 * (verts[maxlam]+verts[(maxlam+1)%3]);
      Vec<2> p2 = { 1.0/3, 1.0/3 };
      Vec<2> p3 = 0.5 * (verts[maxlam]+verts[(maxlam+2)%3]);

      Vec<2> dpdxi = eta*(p2-p3)+(1-eta)*(p1-p0);
      Vec<2> dpdeta = xi*(p2-p1)+(1-xi)*(p3-p0);

      Mat<2,2> dxy_dxieta;
      dxy_dxieta.Col(0) = dpdxi;
      dxy_dxieta.Col(1) = dpdeta;
      Mat<2,2> trafo = 1.0/Det(dxy_dxieta)*dxy_dxieta; // Piola
      
      auto assign_shapex = [&](int nr, int ix, int iy)
        {
          double xi = IR[ix](0);
          double eta = IR[iy](0);
          Vec<2> dpdxi = eta*(p2-p3)+(1-eta)*(p1-p0);
          Vec<2> dpdeta = xi*(p2-p1)+(1-xi)*(p3-p0);

          Mat<2,2> dxy_dxieta;
          dxy_dxieta.Col(0) = dpdxi;
          dxy_dxieta.Col(1) = dpdeta;
          Mat<2,2> trafonode = 1.0/Det(dxy_dxieta)*dxy_dxieta; // Piola
          shape.Row(nr) = trafonode*Vec<2>(polx[ix]*poly[iy], 0);
        };
      auto assign_shapey = [&](int nr, int ix, int iy)
        {
          double xi = IR[ix](0);
          double eta = IR[iy](0);
          Vec<2> dpdxi = eta*(p2-p3)+(1-eta)*(p1-p0);
          Vec<2> dpdeta = xi*(p2-p1)+(1-xi)*(p3-p0);

          Mat<2,2> dxy_dxieta;
          dxy_dxieta.Col(0) = dpdxi;
          dxy_dxieta.Col(1) = dpdeta;
          Mat<2,2> trafonode = 1.0/Det(dxy_dxieta)*dxy_dxieta; // Piola
          shape.Row(nr) = trafonode*Vec<2>(0, -polx[ix]*poly[iy]);
        };

      int ii = 0;

      // x-edge
      ii = (maxlam+2)%3 * ndof_edge;
      for (int i = 0; i < poly.Size(); i++)      
        assign_shapex(ii++, 0, i);
      
      ii = (maxlam+1)%3 * ndof_edge;
      for (int i = 0; i < polx.Size(); i++)
        assign_shapey(ii++, i, 0);
      
      ii = 3*ndof_edge + maxlam*ndof_quad;
      for (int i = 1; i < polx.Size(); i++)
        for (int j = 0; j < poly.Size(); j++)
          assign_shapex(ii++, i, j);
      for (int j = 1; j < poly.Size(); j++)
        for (int i = 0; i < polx.Size(); i++)
          assign_shapey(ii++, i, j);
    }
    void HDivPrimalCellTrig:: CalcCovShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<> shape) const
    {

      //
      // int ndof_edge = order+1;
      // int ndof_quad = 2*sqr(order+1);
      int ndof_edge = IR.Size();
      int ndof_quad = 2*IR.Size()*(IR.Size()-1);
      
      ArrayMem<double, 20> polx(IR.Size()), poly(IR.Size());      
      
      double lam[] = { ip(0), ip(1), 1-ip(0)-ip(1) };
      int maxlam = 0;
      for (int j = 1; j < 3; j++)
        if (lam[j] > lam[maxlam])
          maxlam = j;
      
      //shape = 0;
      shape.AddSize(ndof, 2) = 0;

      double lx = lam[(maxlam+1)%3];
      double ly = lam[(maxlam+2)%3];
      double xi  = (3+2*lx-2*ly)/2 - sqrt( sqr( (2*ly-2*lx-3)/2 ) - 6*lx );      
      double eta = (3+2*ly-2*lx)/2 - sqrt( sqr( (2*lx-2*ly-3)/2 ) - 6*ly );

      LagrangePolynomials(xi, IR, polx);
      LagrangePolynomials(eta, IR, poly);

      Vec<2> verts[3] = { { 1, 0 }, { 0, 1 }, { 0, 0 } };      
      Vec<2> p0 = verts[maxlam];
      Vec<2> p1 = 0.5 * (verts[maxlam]+verts[(maxlam+1)%3]);
      Vec<2> p2 = { 1.0/3, 1.0/3 };
      Vec<2> p3 = 0.5 * (verts[maxlam]+verts[(maxlam+2)%3]);

      Vec<2> dpdxi = eta*(p2-p3)+(1-eta)*(p1-p0);
      Vec<2> dpdeta = xi*(p2-p1)+(1-xi)*(p3-p0);

      Mat<2,2> dxy_dxieta;
      dxy_dxieta.Col(0) = dpdxi;
      dxy_dxieta.Col(1) = dpdeta;
      Mat<2,2> trafo = 1.0/Det(dxy_dxieta)*dxy_dxieta; // Piola
      
      auto assign_shapex = [&](int nr, int ix, int iy)
        {
          double xi = IR[ix](0);
          double eta = IR[iy](0);
          Vec<2> dpdxi = eta*(p2-p3)+(1-eta)*(p1-p0);
          Vec<2> dpdeta = xi*(p2-p1)+(1-xi)*(p3-p0);

          Mat<2,2> dxy_dxieta;
          dxy_dxieta.Col(0) = dpdxi;
          dxy_dxieta.Col(1) = dpdeta;
          //Mat<2,2> trafonode = Trans(Inv(dxy_dxieta)); // Cov 
          Mat<2,2> trafonode = 1.0/Det(dxy_dxieta)*dxy_dxieta; // Piola
          shape.Row(nr) = trafonode*Vec<2>(polx[ix]*poly[iy], 0);
        };
      auto assign_shapey = [&](int nr, int ix, int iy)
        {
          double xi = IR[ix](0);
          double eta = IR[iy](0);
          Vec<2> dpdxi = eta*(p2-p3)+(1-eta)*(p1-p0);
          Vec<2> dpdeta = xi*(p2-p1)+(1-xi)*(p3-p0);

          Mat<2,2> dxy_dxieta;
          dxy_dxieta.Col(0) = dpdxi;
          dxy_dxieta.Col(1) = dpdeta;
          //Mat<2,2> trafonode = Trans(Inv(dxy_dxieta)); // Cov 
          Mat<2,2> trafonode = 1.0/Det(dxy_dxieta)*dxy_dxieta; // Piola
          shape.Row(nr) = trafonode*Vec<2>(0, -polx[ix]*poly[iy]);
        };

      int ii = 0;

      // x-edge
      ii = (maxlam+2)%3 * ndof_edge;
      for (int i = 0; i < poly.Size(); i++)      
        assign_shapex(ii++, 0, i);
      
      ii = (maxlam+1)%3 * ndof_edge;
      for (int i = 0; i < polx.Size(); i++)
        assign_shapey(ii++, i, 0);
      
      ii = 3*ndof_edge + maxlam*ndof_quad;
      for (int i = 1; i < polx.Size(); i++)
        for (int j = 0; j < poly.Size(); j++)
          assign_shapex(ii++, i, j);
      for (int j = 1; j < poly.Size(); j++)
        for (int i = 0; i < polx.Size(); i++)
          assign_shapey(ii++, i, j);
    }

    void HDivPrimalCellTrig:: CalcAltCovShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<> shape) const
    {

      //
      // int ndof_edge = order+1;
      // int ndof_quad = 2*sqr(order+1);
      int ndof_edge = IR.Size();
      int ndof_quad = 2*IR.Size()*(IR.Size()-1);
      
      ArrayMem<double, 20> polx(IR.Size()), poly(IR.Size());      
      
      double lam[] = { ip(0), ip(1), 1-ip(0)-ip(1) };
      int maxlam = 0;
      for (int j = 1; j < 3; j++)
        if (lam[j] > lam[maxlam])
          maxlam = j;
      
      //shape = 0;
      shape.AddSize(ndof, 2) = 0;

      double lx = lam[(maxlam+1)%3];
      double ly = lam[(maxlam+2)%3];
      double xi  = (3+2*lx-2*ly)/2 - sqrt( sqr( (2*ly-2*lx-3)/2 ) - 6*lx );      
      double eta = (3+2*ly-2*lx)/2 - sqrt( sqr( (2*lx-2*ly-3)/2 ) - 6*ly );

      LagrangePolynomials(xi, IR, polx);
      LagrangePolynomials(eta, IR, poly);

      Vec<2> verts[3] = { { 1, 0 }, { 0, 1 }, { 0, 0 } };      
      Vec<2> p0 = verts[maxlam];
      Vec<2> p1 = 0.5 * (verts[maxlam]+verts[(maxlam+1)%3]);
      Vec<2> p2 = { 1.0/3, 1.0/3 };
      Vec<2> p3 = 0.5 * (verts[maxlam]+verts[(maxlam+2)%3]);

      Vec<2> dpdxi = eta*(p2-p3)+(1-eta)*(p1-p0);
      Vec<2> dpdeta = xi*(p2-p1)+(1-xi)*(p3-p0);

      Mat<2,2> dxy_dxieta;
      dxy_dxieta.Col(0) = dpdxi;
      dxy_dxieta.Col(1) = dpdeta;
      Mat<2,2> trafo = 1.0/Det(dxy_dxieta)*dxy_dxieta; // Piola
      
      auto assign_shapex = [&](int nr, int ix, int iy)
        {
          double xi = IR[ix](0);
          double eta = IR[iy](0);
          Vec<2> dpdxi = eta*(p2-p3)+(1-eta)*(p1-p0);
          Vec<2> dpdeta = xi*(p2-p1)+(1-xi)*(p3-p0);

          Mat<2,2> dxy_dxieta;
          dxy_dxieta.Col(0) = dpdxi;
          dxy_dxieta.Col(1) = dpdeta;
          Mat<2,2> trafonode = Trans(Inv(dxy_dxieta)); // Cov 
          shape.Row(nr) = trafonode*Vec<2>(polx[ix]*poly[iy], 0);
        };
      auto assign_shapey = [&](int nr, int ix, int iy)
        {
          double xi = IR[ix](0);
          double eta = IR[iy](0);
          Vec<2> dpdxi = eta*(p2-p3)+(1-eta)*(p1-p0);
          Vec<2> dpdeta = xi*(p2-p1)+(1-xi)*(p3-p0);

          Mat<2,2> dxy_dxieta;
          dxy_dxieta.Col(0) = dpdxi;
          dxy_dxieta.Col(1) = dpdeta;
          Mat<2,2> trafonode = Trans(Inv(dxy_dxieta)); // Cov 
          shape.Row(nr) = trafonode*Vec<2>(0, -polx[ix]*poly[iy]);
        };

      int ii = 0;

      // x-edge
      ii = (maxlam+2)%3 * ndof_edge;
      for (int i = 0; i < poly.Size(); i++)      
        assign_shapex(ii++, 0, i);
      
      ii = (maxlam+1)%3 * ndof_edge;
      for (int i = 0; i < polx.Size(); i++)
        assign_shapey(ii++, i, 0);
      
      ii = 3*ndof_edge + maxlam*ndof_quad;
      for (int i = 1; i < polx.Size(); i++)
        for (int j = 0; j < poly.Size(); j++)
          assign_shapex(ii++, i, j);
      for (int j = 1; j < poly.Size(); j++)
        for (int i = 0; i < polx.Size(); i++)
          assign_shapey(ii++, i, j);
    }


    void HDivPrimalCellTrig::
    CalcDivShape (const IntegrationPoint & ip, 
                               SliceVector<> divshape) const
    {
      // int ndof_edge = order+1;
      // int ndof_quad = 2*sqr(order+1);
      int ndof_edge = IR.Size();
      int ndof_quad = 2*IR.Size()*(IR.Size()-1);
      
      ArrayMem<AutoDiff<1>, 20> Dpolx(IR.Size()), Dpoly(IR.Size());
      ArrayMem<double, 20> polx(IR.Size()), poly(IR.Size());      
      
      double lam[] = { ip(0), ip(1), 1-ip(0)-ip(1) };
      int maxlam = 0;
      for (int j = 1; j < 3; j++)
        if (lam[j] > lam[maxlam])
          maxlam = j;

      divshape = 0;
      //divshape.AddSize(ndof, 2) = 0;

      double lx = lam[(maxlam+1)%3];
      double ly = lam[(maxlam+2)%3];
      double xi  = (3+2*lx-2*ly)/2 - sqrt( sqr( (2*ly-2*lx-3)/2 ) - 6*lx );      
      double eta = (3+2*ly-2*lx)/2 - sqrt( sqr( (2*lx-2*ly-3)/2 ) - 6*ly );

      AutoDiff<1> adxi(xi, 0);
      AutoDiff<1> adeta(eta, 0);
      
      LagrangePolynomials(adxi, IR, Dpolx);
      LagrangePolynomials(adeta, IR, Dpoly);
      LagrangePolynomials(xi, IR, polx);
      LagrangePolynomials(eta, IR, poly);

      Vec<2> verts[3] = { { 1, 0 }, { 0, 1 }, { 0, 0 } };      
      Vec<2> p0 = verts[maxlam];
      Vec<2> p1 = 0.5 * (verts[maxlam]+verts[(maxlam+1)%3]);
      Vec<2> p2 = { 1.0/3, 1.0/3 };
      Vec<2> p3 = 0.5 * (verts[maxlam]+verts[(maxlam+2)%3]);

      Vec<2> dpdxi = eta*(p2-p3)+(1-eta)*(p1-p0);
      Vec<2> dpdeta = xi*(p2-p1)+(1-xi)*(p3-p0);

      Mat<2,2> dxy_dxieta;
      dxy_dxieta.Col(0) = dpdxi;
      dxy_dxieta.Col(1) = dpdeta;
      double trafo = 1.0/Det(dxy_dxieta);
      
      auto assign_shapex = [trafo, divshape, &Dpolx, &poly](int nr, int ix, int iy)
        {
          divshape(nr) = trafo*Dpolx[ix].DValue(0)*poly[iy];
        };
      auto assign_shapey = [trafo, divshape, &polx, &Dpoly](int nr, int ix, int iy)
        {
          divshape(nr) = trafo*(-polx[ix]*Dpoly[iy].DValue(0));
        };

      int ii = 0;

      // x-edge
      ii = (maxlam+2)%3 * ndof_edge;
      for (int i = 0; i < poly.Size(); i++)              
        assign_shapex(ii++, 0, i);
      
      ii = (maxlam+1)%3 * ndof_edge;
      for (int i = 0; i < polx.Size(); i++)      
        assign_shapey(ii++, i, 0);
      
      ii = 3*ndof_edge + maxlam*ndof_quad;
      for (int i = 1; i < Dpolx.Size(); i++)
        for (int j = 0; j < Dpoly.Size(); j++)
          assign_shapex(ii++, i, j);
      for (int j = 1; j < Dpoly.Size(); j++)
        for (int i = 0; i < polx.Size(); i++)
          assign_shapey(ii++, i, j);
    }
  
    void HDivPrimalCellTet:: CalcShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<> shape) const
    {
      double lam[] = { ip(0), ip(1), ip(2), 1-ip(0)-ip(1)-ip(2) };
      int maxlam = PosMax(lam);

      //shape = 0;
      shape.AddSize(ndof, 3) = 0;

      int minvi = (maxlam+1)%4;
      int maxvi = minvi;
      for (int i = 0; i < 4; i++)
        if (i != maxlam)
          {
            if (vnums[i] < vnums[minvi]) minvi = i;
            if (vnums[i] > vnums[maxvi]) maxvi = i;
          }
      int midvi = 6-maxlam-minvi-maxvi;

      int vdir[4];
      vdir[maxlam] = -1;
      vdir[minvi] = 0;
      vdir[midvi] = 1;
      vdir[maxvi] = 2;
      
      Vec<3> x(lam[minvi], lam[midvi], lam[maxvi]);
      Vec<3> xi = MapTet2Hex (x);

      Mat<3,3> F = DMapHex2Tet(xi);
      Mat<3,3> F2;    // trafo from vertex permutation
      Vec<3> verts[] = { Vec<3>(1,0,0), Vec<3>(0,1,0), Vec<3> (0,0,1), Vec<3> (0,0,0) };
      F2.Col(0) = verts[minvi]-verts[maxlam];
      F2.Col(1) = verts[midvi]-verts[maxlam];
      F2.Col(2) = verts[maxvi]-verts[maxlam];
      
      Mat<3> trafo = 1.0/fabs(Det(F2*F))*(F2*F);

      int nd = IR.Size();

      
      ArrayMem<double, 20> polxi(nd), poleta(nd), polzeta(nd);
      LagrangePolynomials(xi(0), IR, polxi);
      LagrangePolynomials(xi(1), IR, poleta);
      LagrangePolynomials(xi(2), IR, polzeta);

      auto assign =  [&](int nr, IVec<3> i, int dir,bool flipsign = false)
        {
          double sig = 1;
          if (flipsign)
            sig = -1;
          switch (dir)
            {
            case 0:
            {
              shape.Row(nr) = sig*trafo*Vec<3>(polxi[i[0]]*poleta[i[1]]*polzeta[i[2]], 0, 0);
              break;
            }
            case 1:
            {
              shape.Row(nr) = sig*trafo*Vec<3>(0, polxi[i[0]]*poleta[i[1]]*polzeta[i[2]], 0);
              break;
            }
            case 2:
            {
              shape.Row(nr) = sig*trafo*Vec<3>(0, 0, polxi[i[0]]*poleta[i[1]]*polzeta[i[2]]);
              break;
            }
            break;
            }
        };

      int ii = 0;
      

      // internal faces normal to direction of vertices v1, v2
      for (int v1 = 0; v1 < 4; v1++)
        for (int v2 = 0; v2 < v1; v2++)
          {
            if (v1 != maxlam && v2 != maxlam)
              {
                IVec<3> ind = { 0, 0, 0 };
                int dirv1 = vdir[v1];
                int dirv2 = vdir[v2];
                int kk = ii;
                for (int l = 0; l < nd; l++)
                  for (int k = 0; k < nd; k++)
                    {
                      ind[dirv1] = l; 
                      ind[dirv2] = k;
                      assign(kk++, ind, 3-dirv1-dirv2,6-v1-v2-maxlam<maxlam);
                      //assign(kk++, ind, 3-dirv1-dirv2);
                    }
              }
            ii += nd*nd;
          }
      ii += maxlam*3*sqr(nd)*(nd-1);


      // remaining dofs
      for (int i = 0; i < nd; i++)
        for (int j = 0; j < nd; j++)
          for (int k = 1; k < nd; k++)
            {
              assign(ii++, { k, i, j }, 0);
              assign(ii++, { i, k, j }, 1);
              assign(ii++, { i, j, k }, 2);
            }

    }


    void HDivPrimalCellTet:: CalcAltShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<> shape) const
    {
      double lam[] = { ip(0), ip(1), ip(2), 1-ip(0)-ip(1)-ip(2) };
      int maxlam = PosMax(lam);

      //shape = 0;
      shape.AddSize(ndof, 3) = 0;

      int minvi = (maxlam+1)%4;
      int maxvi = minvi;
      for (int i = 0; i < 4; i++)
        if (i != maxlam)
          {
            if (vnums[i] < vnums[minvi]) minvi = i;
            if (vnums[i] > vnums[maxvi]) maxvi = i;
          }
      int midvi = 6-maxlam-minvi-maxvi;



      int vdir[4];
      vdir[maxlam] = -1;
      vdir[minvi] = 0;
      vdir[midvi] = 1;
      vdir[maxvi] = 2;

      //cout << "global numbers " << vnums[0] << vnums[1] << vnums[2] << vnums[3] << endl;
      //cout << "maxlam " << maxlam << "  minvi, midvi,maxvi" << minvi<< midvi << maxvi << endl;
      //cout << "vdir " << vdir[0] << vdir[1] << vdir[2] << vdir[3] << endl;

      
      Vec<3> x(lam[minvi], lam[midvi], lam[maxvi]);
      Vec<3> xi = MapTet2Hex (x);

      Mat<3,3> F = DMapHex2Tet(xi);
      Mat<3,3> F2;    // trafo from vertex permutation
      Vec<3> verts[] = { Vec<3>(1,0,0), Vec<3>(0,1,0), Vec<3> (0,0,1), Vec<3> (0,0,0) };
      F2.Col(0) = verts[minvi]-verts[maxlam];
      F2.Col(1) = verts[midvi]-verts[maxlam];
      F2.Col(2) = verts[maxvi]-verts[maxlam];
      
      Mat<3> trafo = 1.0/fabs(Det(F2*F))*(F2*F);

      int nd = IR.Size();

      
      ArrayMem<double, 20> polxi(nd), poleta(nd), polzeta(nd);
      LagrangePolynomials(xi(0), IR, polxi);
      LagrangePolynomials(xi(1), IR, poleta);
      LagrangePolynomials(xi(2), IR, polzeta);


      auto assign =  [&](int nr, IVec<3> i, int dir,bool flipsign = false)
        {
          //cout << "assigning nr " << nr << " i= " << i << " dir = " << dir << " flip " << flipsign << endl;
          Vec<3> xinode(IR[i[0]](0), IR[i[1]](0), IR[i[2]](0));
          Mat<3,3> Fnode = DMapHex2Tet(xinode);          
          Mat<3> trafonode = 1.0/fabs(Det(F2*Fnode))*(F2*Fnode);

          //cout << "xinode " << xinode[0] << xinode[1] << xinode[2] <<  endl;
          
          if (flipsign)
            trafonode *=-1.;
          switch (dir)
            {
            case 0:
            {
              shape.Row(nr) = trafonode*Vec<3>(polxi[i[0]]*poleta[i[1]]*polzeta[i[2]], 0, 0);
              break;
            }
            case 1:
            {
              shape.Row(nr) = trafonode*Vec<3>(0, polxi[i[0]]*poleta[i[1]]*polzeta[i[2]], 0);
              break;
            }
            case 2:
            {
              shape.Row(nr) = trafonode*Vec<3>(0, 0, polxi[i[0]]*poleta[i[1]]*polzeta[i[2]]);
              break;
            }
            //break;
            }
        };

      int ii = 0;
      

      // internal faces normal to direction of vertices v1, v2
      //cout << "maxlam = " << maxlam << endl;
      //cout << "vdir = " << vdir[0] << vdir[1]<< vdir[2] << vdir[3] << endl;
      for (int v1 = 0; v1 < 4; v1++)
        for (int v2 = 0; v2 < v1; v2++)
          {
            if (v1 != maxlam && v2 != maxlam)
              {
                //cout << "v1 = " << v1 << "v2 = " << v2 << endl;
              
                IVec<3> ind = { 0, 0, 0 };
                int dirv1 = vdir[v1];
                int dirv2 = vdir[v2];
                int kk = ii;
                for (int l = 0; l < nd; l++)
                  for (int k = 0; k < nd; k++)
                    {
                      ind[dirv1] = l; 
                      ind[dirv2] = k;
                      assign(kk++, ind, 3-dirv1-dirv2,maxlam>6-v1-v2-maxlam);
                      //assign(kk++, ind, 3-dirv1-dirv2);
                    }
              }
            ii += nd*nd;
          }
      ii += maxlam*3*sqr(nd)*(nd-1);


      // remaining dofs
      for (int i = 0; i < nd; i++)
        for (int j = 0; j < nd; j++)
          for (int k = 1; k < nd; k++)
            {
              assign(ii++, { k, i, j }, 0);
              assign(ii++, { i, k, j }, 1);
              assign(ii++, { i, j, k }, 2);
            }

    }
    void HDivPrimalCellTet:: CalcCovShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<> shape) const
    {
      double lam[] = { ip(0), ip(1), ip(2), 1-ip(0)-ip(1)-ip(2) };
      int maxlam = PosMax(lam);

      //shape = 0;
      shape.AddSize(ndof, 3) = 0;

      int minvi = (maxlam+1)%4;
      int maxvi = minvi;
      for (int i = 0; i < 4; i++)
        if (i != maxlam)
          {
            if (vnums[i] < vnums[minvi]) minvi = i;
            if (vnums[i] > vnums[maxvi]) maxvi = i;
          }
      int midvi = 6-maxlam-minvi-maxvi;



      int vdir[4];
      vdir[maxlam] = -1;
      vdir[minvi] = 0;
      vdir[midvi] = 1;
      vdir[maxvi] = 2;

      //cout << "global numbers " << vnums[0] << vnums[1] << vnums[2] << vnums[3] << endl;
      //cout << "maxlam " << maxlam << "  minvi, midvi,maxvi" << minvi<< midvi << maxvi << endl;
      //cout << "vdir " << vdir[0] << vdir[1] << vdir[2] << vdir[3] << endl;

      
      Vec<3> x(lam[minvi], lam[midvi], lam[maxvi]);
      Vec<3> xi = MapTet2Hex (x);

      Mat<3,3> F = DMapHex2Tet(xi);
      Mat<3,3> F2;    // trafo from vertex permutation
      Vec<3> verts[] = { Vec<3>(1,0,0), Vec<3>(0,1,0), Vec<3> (0,0,1), Vec<3> (0,0,0) };
      F2.Col(0) = verts[minvi]-verts[maxlam];
      F2.Col(1) = verts[midvi]-verts[maxlam];
      F2.Col(2) = verts[maxvi]-verts[maxlam];
      
      Mat<3> trafo = 1.0/fabs(Det(F2*F))*(F2*F);

      int nd = IR.Size();

      
      ArrayMem<double, 20> polxi(nd), poleta(nd), polzeta(nd);
      LagrangePolynomials(xi(0), IR, polxi);
      LagrangePolynomials(xi(1), IR, poleta);
      LagrangePolynomials(xi(2), IR, polzeta);

      auto assign =  [&](int nr, IVec<3> i, int dir,bool flipsign = false)
        {
          //cout << "assigning nr " << nr << " i= " << i << " dir = " << dir << " flip " << flipsign << endl;
          Vec<3> xinode(IR[i[0]](0), IR[i[1]](0), IR[i[2]](0));
          Mat<3,3> Fnode = DMapHex2Tet(xinode);          
          Mat<3> trafonode = 1.0/fabs(Det(F2*Fnode))*(F2*Fnode);

          //cout << "xinode " << xinode[0] << xinode[1] << xinode[2] <<  endl;

          if (flipsign)
            trafonode*=-1.;
          switch (dir)
            {
            case 0:
            {
              shape.Row(nr) = trafonode*Vec<3>(polxi[i[0]]*poleta[i[1]]*polzeta[i[2]], 0, 0);
              break;
            }
            case 1:
            {
              shape.Row(nr) = trafonode*Vec<3>(0, polxi[i[0]]*poleta[i[1]]*polzeta[i[2]], 0);
              break;
            }
            case 2:
            {
              shape.Row(nr) = trafonode*Vec<3>(0, 0, polxi[i[0]]*poleta[i[1]]*polzeta[i[2]]);
              break;
            }
            //break;
            }
        };

      int ii = 0;
      

      // internal faces normal to direction of vertices v1, v2
      //cout << "maxlam = " << maxlam << endl;
      //cout << "vdir = " << vdir[0] << vdir[1]<< vdir[2] << vdir[3] << endl;
      for (int v1 = 0; v1 < 4; v1++)
        for (int v2 = 0; v2 < v1; v2++)
          {
            if (v1 != maxlam && v2 != maxlam)
              {
                //cout << "v1 = " << v1 << "v2 = " << v2 << endl;
              
                IVec<3> ind = { 0, 0, 0 };
                int dirv1 = vdir[v1];
                int dirv2 = vdir[v2];
                int kk = ii;
                for (int l = 0; l < nd; l++)
                  for (int k = 0; k < nd; k++)
                    {
                      ind[dirv1] = l; 
                      ind[dirv2] = k;
                      assign(kk++, ind, 3-dirv1-dirv2,maxlam>6-v1-v2-maxlam);
                      //assign(kk++, ind, 3-dirv1-dirv2);
                    }
              }
            ii += nd*nd;
          }
      ii += maxlam*3*sqr(nd)*(nd-1);


      // remaining dofs
      for (int i = 0; i < nd; i++)
        for (int j = 0; j < nd; j++)
          for (int k = 1; k < nd; k++)
            {
              assign(ii++, { k, i, j }, 0);
              assign(ii++, { i, k, j }, 1);
              assign(ii++, { i, j, k }, 2);
            }

    }

    void HDivPrimalCellTet::
    CalcDivShape (const IntegrationPoint & ip, 
                               SliceVector<> divshape) const
    {
      double lam[] = { ip(0), ip(1), ip(2), 1-ip(0)-ip(1)-ip(2) };
      int maxlam = PosMax(lam);

      divshape = 0;
      //divshape.AddSize(ndof, 3) = 0;

      int minvi = (maxlam+1)%4;
      int maxvi = minvi;
      for (int i = 0; i < 4; i++)
        if (i != maxlam)
          {
            if (vnums[i] < vnums[minvi]) minvi = i;
            if (vnums[i] > vnums[maxvi]) maxvi = i;
          }
      int midvi = 6-maxlam-minvi-maxvi;

      int vdir[4];
      vdir[maxlam] = -1;
      vdir[minvi] = 0;
      vdir[midvi] = 1;
      vdir[maxvi] = 2;
      
      Vec<3> x(lam[minvi], lam[midvi], lam[maxvi]);
      Vec<3> xi = MapTet2Hex (x);

      Mat<3,3> F = DMapHex2Tet(xi);
      Mat<3,3> F2;    // trafo from vertex permutation
      Vec<3> verts[] = { Vec<3>(1,0,0), Vec<3>(0,1,0), Vec<3> (0,0,1), Vec<3> (0,0,0) };
      F2.Col(0) = verts[minvi]-verts[maxlam];
      F2.Col(1) = verts[midvi]-verts[maxlam];
      F2.Col(2) = verts[maxvi]-verts[maxlam];
      
      double trafo = 1.0/fabs(Det(F2*F));

      int nd = IR.Size();

      ArrayMem<AutoDiff<1>, 20> Dpolxi(nd), Dpoleta(nd), Dpolzeta(nd);
      ArrayMem<double, 20> polxi(nd), poleta(nd), polzeta(nd);         
      LagrangePolynomials(AutoDiff<1>(xi(0),0), IR, Dpolxi);
      LagrangePolynomials(AutoDiff<1>(xi(1),0), IR, Dpoleta);
      LagrangePolynomials(AutoDiff<1>(xi(2),0), IR, Dpolzeta);

      LagrangePolynomials(xi(0), IR, polxi);
      LagrangePolynomials(xi(1), IR, poleta);
      LagrangePolynomials(xi(2), IR, polzeta);

      auto assign =  [&](int nr, IVec<3> i, int dir,bool flipsign = false)
        {
          switch (dir)
            {
            case 0:
            {
              divshape(nr) = trafo*Dpolxi[i[0]].DValue(0)*poleta[i[1]]*polzeta[i[2]];
              break;
            }
            case 1:
            {
              divshape(nr) = trafo*polxi[i[0]]*Dpoleta[i[1]].DValue(0)*polzeta[i[2]];
              break;
            }
            case 2:
            {
              divshape(nr) = trafo*polxi[i[0]]*poleta[i[1]]*Dpolzeta[i[2]].DValue(0);
              break;
            }
            break;
            }
          if (flipsign)
            divshape(nr)*=-1.;
        };

      int ii = 0;
      

      // internal faces normal to direction of vertices v1, v2
      for (int v1 = 0; v1 < 4; v1++)
        for (int v2 = 0; v2 < v1; v2++)
          {
            if (v1 != maxlam && v2 != maxlam)
              {
                IVec<3> ind = { 0, 0, 0 };
                int dirv1 = vdir[v1];
                int dirv2 = vdir[v2];
                int kk = ii;
                for (int l = 0; l < nd; l++)
                  for (int k = 0; k < nd; k++)
                    {
                      ind[dirv1] = l; 
                      ind[dirv2] = k;
                      assign(kk++, ind, 3-dirv1-dirv2,6-v1-v2-maxlam<maxlam);
                      //assign(kk++, ind, 3-dirv1-dirv2);
                    }
              }
            ii += nd*nd;
          }
      ii += maxlam*3*sqr(nd)*(nd-1);


      // remaining dofs
      for (int i = 0; i < nd; i++)
        for (int j = 0; j < nd; j++)
          for (int k = 1; k < nd; k++)
            {
              assign(ii++, { k, i, j }, 0);
              assign(ii++, { i, k, j }, 1);
              assign(ii++, { i, j, k }, 2);
            }

    }
  
    void HDivPrimalCellTet::
    CalcAltDivShape (const IntegrationPoint & ip, 
                               BareSliceMatrix<> divshape) const
    {
      double lam[] = { ip(0), ip(1), ip(2), 1-ip(0)-ip(1)-ip(2) };
      int maxlam = PosMax(lam);

      divshape.AddSize(ndof, 1) = 0;
      //divshape = 0;

      int minvi = (maxlam+1)%4;
      int maxvi = minvi;
      for (int i = 0; i < 4; i++)
        if (i != maxlam)
          {
            if (vnums[i] < vnums[minvi]) minvi = i;
            if (vnums[i] > vnums[maxvi]) maxvi = i;
          }
      int midvi = 6-maxlam-minvi-maxvi;

      int vdir[4];
      vdir[maxlam] = -1;
      vdir[minvi] = 0;
      vdir[midvi] = 1;
      vdir[maxvi] = 2;
      
      Vec<3> x(lam[minvi], lam[midvi], lam[maxvi]);
      Vec<3> xi = MapTet2Hex (x);

      Mat<3,3> F = DMapHex2Tet(xi);
      Mat<3,3> F2;    // trafo from vertex permutation
      Vec<3> verts[] = { Vec<3>(1,0,0), Vec<3>(0,1,0), Vec<3> (0,0,1), Vec<3> (0,0,0) };
      F2.Col(0) = verts[minvi]-verts[maxlam];
      F2.Col(1) = verts[midvi]-verts[maxlam];
      F2.Col(2) = verts[maxvi]-verts[maxlam];
      
      Mat<3,3> trafo = Inv(F2*F);

      int nd = IR.Size();

      ArrayMem<AutoDiff<1>, 20> Dpolxi(nd), Dpoleta(nd), Dpolzeta(nd);         
      LagrangePolynomials(AutoDiff<1>(xi(0),0), IR, Dpolxi);
      LagrangePolynomials(AutoDiff<1>(xi(1),0), IR, Dpoleta);
      LagrangePolynomials(AutoDiff<1>(xi(2),0), IR, Dpolzeta);

      auto assign =  [&](int nr, IVec<3> i, int dir,bool flipsign = false)
        {
          Vec<3> xinode(IR[i[0]](0), IR[i[1]](0), IR[i[2]](0));
          Mat<3,3> Fnode = DMapHex2Tet(xinode);          
          Mat<3,3> trafonode = 1.0/fabs(Det(F2*Fnode))*(F2*Fnode);
          //trafonode = 1.0/fabs(Det(F2*F))*(F2*F); //correct piola

          Mat<3,3> Dshapeorig = 0.;
          switch (dir)
            {
            case 0:
            {
            Dshapeorig.Col(0) = Trans(trafo)*Vec<3>(Dpolxi[i[0]].DValue(0)*Dpoleta[i[1]].Value()*Dpolzeta[i[2]].Value(),
              Dpolxi[i[0]].Value()*Dpoleta[i[1]].DValue(0)*Dpolzeta[i[2]].Value(),
              Dpolxi[i[0]].Value()*Dpoleta[i[1]].Value()*Dpolzeta[i[2]].DValue(0));
              break;
            }
            case 1:
            {
            Dshapeorig.Col(1) = Trans(trafo)*Vec<3>(Dpolxi[i[0]].DValue(0)*Dpoleta[i[1]].Value()*Dpolzeta[i[2]].Value(),
              Dpolxi[i[0]].Value()*Dpoleta[i[1]].DValue(0)*Dpolzeta[i[2]].Value(),
              Dpolxi[i[0]].Value()*Dpoleta[i[1]].Value()*Dpolzeta[i[2]].DValue(0));
              break;
            }
            case 2:
            {
            Dshapeorig.Col(2) = Trans(trafo)*Vec<3>(Dpolxi[i[0]].DValue(0)*Dpoleta[i[1]].Value()*Dpolzeta[i[2]].Value(),
              Dpolxi[i[0]].Value()*Dpoleta[i[1]].DValue(0)*Dpolzeta[i[2]].Value(),
              Dpolxi[i[0]].Value()*Dpoleta[i[1]].Value()*Dpolzeta[i[2]].DValue(0));
              break;
            }
            }
          if (flipsign)
            Dshapeorig*=-1.;
          Mat<3,3> Dshape = trafonode*Trans(Dshapeorig);
          divshape(nr) = Trace(Dshape);

        };

      int ii = 0;
      

      // internal faces normal to direction of vertices v1, v2
      for (int v1 = 0; v1 < 4; v1++)
        for (int v2 = 0; v2 < v1; v2++)
          {
            if (v1 != maxlam && v2 != maxlam)
              {
                IVec<3> ind = { 0, 0, 0 };
                int dirv1 = vdir[v1];
                int dirv2 = vdir[v2];
                int kk = ii;
                for (int l = 0; l < nd; l++)
                  for (int k = 0; k < nd; k++)
                    {
                      ind[dirv1] = l; 
                      ind[dirv2] = k;
                      assign(kk++, ind, 3-dirv1-dirv2,6-v1-v2-maxlam<maxlam);
                      //assign(kk++, ind, 3-dirv1-dirv2);
                    }
              }
            ii += nd*nd;
          }
      ii += maxlam*3*sqr(nd)*(nd-1);


      // remaining dofs
      for (int i = 0; i < nd; i++)
        for (int j = 0; j < nd; j++)
          for (int k = 1; k < nd; k++)
            {
              assign(ii++, { k, i, j }, 0);
              assign(ii++, { i, k, j }, 1);
              assign(ii++, { i, j, k }, 2);
            }

    }


  HDivPrimalCells::
    HDivPrimalCells (shared_ptr<MeshAccess> ama, const Flags & flags)
      : FESpace (ama, flags)
    {

      //DefineDefineFlag("uniform_order");
      if (ma->GetDimension()==2)
      {
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHDiv<2>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpDivHDiv<2>>>();
        additional_evaluators.Set ("altshape", make_shared<T_DifferentialOperator<DiffOpAltShapeHDiv<2>>> ());
      }
      if (ma->GetDimension()==3)
      {
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHDiv<3>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpDivHDiv<3>>>();
        additional_evaluators.Set ("altshape", make_shared<T_DifferentialOperator<DiffOpAltShapeHDiv<3>>> ());
      }

      IR = IntegrationRule();
      Array<double> xi, wi;
      ComputeGaussRadauRule (order+1, xi, wi);
      for (auto i : Range(xi))
        IR.Append (IntegrationPoint (1-xi[i]-1e-12, 0, 0, wi[i]));
    }


    // organize the FESpace, called after every mesh update
    void HDivPrimalCells::Update() 
    {
      FESpace::Update();
      size_t ndof = 0;
      
      first_element_dofs.SetSize(ma->GetNE()+1);
      first_element_dofs[0] = ndof;
      switch (ma->GetDimension())
      {
        case 2:
        {
          for (auto i : Range(ma->GetNE()))
            first_element_dofs[i+1] = ndof += 3*IR.Size()*(2*IR.Size()-1);
          break;
        }
        case 3:
        {
          for (auto i : Range(ma->GetNE()))
            first_element_dofs[i+1] = ndof += 6*IR.Size()*IR.Size()*(2*IR.Size()-1);
          break;
        }
        default:
        {
          throw Exception("only dimensions 2 and 3 implemented");
          break;
        }
      }
      SetNDof(ndof);
    }
    
    
    FiniteElement & HDivPrimalCells::
      GetFE (ElementId ei, Allocator & alloc) const 
    {
      auto ngel = ma->GetElement (ei);
      switch (ngel.GetType())
        {
        case ET_TRIG:
          {
            auto trig = new (alloc) HDivPrimalCellTrig(IR);
            trig->SetVertexNumbers (ngel.vertices);
            return *trig;
          }
        case ET_TET:
          {
            auto tet = new (alloc) HDivPrimalCellTet(IR);
            tet->SetVertexNumbers (ngel.vertices);
            return *tet;
          }
        default:
          throw Exception("element not implemented");
        }
    }
    shared_ptr<BaseMatrix> HDivPrimalCells::
      GetMassOperator(shared_ptr<CoefficientFunction> rho, shared_ptr<Region> defon, LocalHeap & lh) const
    {

      if (ma->GetDimension()==3)
      {
      //if (!uniform_order)
      //  throw Exception("3d not implemented");
      
   static Timer t("GetMassOperator"); RegionTimer reg(t);
    static Timer tint("integrate");
    static Timer tsp("MakeSparse");
      
    HeapReset hr(lh);
    auto irs = this->GetIntegrationRules();
    IntegrationRule ir = std::move(irs[ET_TET]);

    auto & felref = dynamic_cast<const HDivPrimalCellTet&> (GetFE(ElementId(VOL,0), lh));
    Matrix shapes(felref.GetNDof(), 3*ir.Size());
    Matrix shapes_trans(3*ir.Size(), felref.GetNDof());
    
    Array<int> rowind, colind;
    Array<double> values;
    
    Array<short> classnr(ma->GetNE(VOL));
    ma->IterateElements
      (VOL, lh, [&] (auto el, LocalHeap & llh)
       {
         classnr[el.Nr()] = 
           SwitchET<ET_SEGM, ET_TRIG,ET_TET>
           (el.GetType(),
            [el] (auto et) { return ET_trait<et.ElementType()>::GetClassNr(el.Vertices()); });
       });

    TableCreator<size_t> creator;
    for ( ; !creator.Done(); creator++)
      for (auto i : Range(classnr))
        creator.Add (classnr[i], i);
    Table<size_t> table = creator.MoveTable();

    
    for (auto elclass_inds : table)
      {
        if (elclass_inds.Size() == 0) continue;

        HeapReset hr(lh);
        

        size_t nr = elclass_inds[0];
        auto & trafo = ma->GetTrafo(ElementId(VOL,nr), lh);
        auto & felref = dynamic_cast<const HDivPrimalCellTet&> (GetFE(ElementId(VOL,nr), lh));
        
        for (int i = 0; i < ir.Size(); i++)
          felref.CalcShape (ir[i], shapes.Cols(3*i, 3*i+3));
        for (int i = 0; i < shapes.Height(); i++)
          for (int j = 0; j < shapes.Width(); j++)
            if (fabs(shapes(i,j)) < 1e-8)
              shapes(i,j) = 0;
        SmallSparseMatrix spshapes(shapes, 1e-8);
        shapes_trans = Trans(shapes);
        
        tint.Start();

        MyMutex COOmutex;
        ParallelForRange (elclass_inds.Size(), [&](auto myrange)
          {
            auto mylh = lh.Split();
            Array<int> myrowind, mycolind;
            Array<double> myvalues;
            Array<DofId> dofs(felref.GetNDof());
            Matrix rhoi_shapes_trans(3*ir.Size(), felref.GetNDof());
            Matrix<> elmat(felref.GetNDof());
            
            for (auto i : myrange)
              {
                HeapReset hr(mylh);
                auto nr = elclass_inds[i];

                auto & trafo = ma->GetTrafo(ElementId(VOL,nr), mylh);            
                MappedIntegrationRule<3,3> mir(ir, trafo, mylh);
                GetDofNrs(ElementId(VOL,nr), dofs);
                
                for (size_t i = 0; i < mir.Size(); i++)
                  {
                    Mat<3,3> rhoi = Id<3>();
                    Mat<3,3> F = mir[i].GetJacobian();
                    rhoi = Trans(F) * rhoi * F;
                    rhoi *= ir[i].Weight() / mir[i].GetJacobiDet();
                    rhoi_shapes_trans.Rows(3*i, 3*i+3) = rhoi * shapes_trans.Rows(3*i, 3*i+3);
                  }
                
                spshapes.Mult (rhoi_shapes_trans, elmat);
                
                for (int i = 0; i < dofs.Size(); i++)
                  for (int j = 0; j < dofs.Size(); j++)
                    if (fabs(elmat(i,j)) > 1e-10)
                      {
                        myrowind.Append(dofs[i]);
                        mycolind.Append(dofs[j]);
                        myvalues.Append(elmat(i,j));
                      }
              }

            COOmutex.lock();
            rowind += myrowind;
            colind += mycolind;
            values += myvalues;
            COOmutex.unlock();
          });

        tint.Stop();
      }    

    
    
      tsp.Start();
      auto spmat = SparseMatrix<double>::CreateFromCOO(rowind, colind, values, GetNDof(), GetNDof());
      tsp.Stop();
      
      return make_shared<MySuperSparseMatrix> (std::move(*spmat));
    }
    if (ma->GetDimension()==2)
    {

      /*if (!uniform_order)
        return make_shared<ApplyMassHDivPrimalCells>
          (dynamic_pointer_cast<FESpace>(const_cast<HDivPrimalCells*>(this)->shared_from_this()),
           rho,false, defon, lh);
           */

      static Timer t("GetMassOperator"); RegionTimer reg(t);
      static Timer tint("integrate");
      static Timer tsp("MakeSparse");
      
      // assemble a (very) sparse mass matrix
      HeapReset hr(lh);
      auto irs = GetIntegrationRules();
      IntegrationRule ir = std::move(irs[ET_TRIG]);

      auto & felref = dynamic_cast<const HDivPrimalCellTrig&> (GetFE(ElementId(VOL,0), lh));
      Matrix shapes(felref.GetNDof(), 2*ir.Size());
      Matrix rhoi_shapes(felref.GetNDof(), 2*ir.Size());

      for (int i = 0; i < ir.Size(); i++)
        felref.CalcShape (ir[i], shapes.Cols(2*i, 2*i+2));

      for (int i = 0; i < shapes.Height(); i++)
        for (int j = 0; j < shapes.Width(); j++)
          if (fabs(shapes(i,j)) < 1e-8)
            shapes(i,j) = 0;
      rhoi_shapes = shapes;
      
      Array<int> rowind, colind;
      Array<double> values;
      Matrix<> elmat(felref.GetNDof());
      Array<DofId> dofs(felref.GetNDof());

      tint.Start();
      for (size_t nr : Range(ma->GetNE()))
        {
          if (defon && !defon->Mask()[ma->GetElIndex(ElementId(VOL,nr))]) continue;
          
          HeapReset hr(lh);
          MappedIntegrationRule<2,2> mir(ir, ma->GetTrafo(ElementId(VOL,nr), lh), lh);
          GetDofNrs(ElementId(VOL,nr), dofs);

          int dimrho = 1;
          if (rho)
            dimrho = rho->Dimension();
          FlatMatrix<> rhovals(mir.Size(), dimrho, lh);
          if (rho)
            rho->Evaluate (mir, rhovals);
          
          for (size_t i = 0; i < mir.Size(); i++)
            {
              Mat<2,2> rhoi = Id<2>();
              if (rho)
                {
                  if (rho->Dimension() == 1)
                    rhoi = rhovals(i,0) * Id<2>();
                  else
                    {
                      auto row = rhovals.Row(i);
                      rhoi = row.AsMatrix(2,2);
                      // rhoi = rhovals.Row(i).AsMatrix(2,2);   // needs ngsolve>=July30
                    }
                }

              
              Mat<2,2> F = mir[i].GetJacobian();
              rhoi = Trans(F) * rhoi * F;
              rhoi *= ir[i].Weight() / mir[i].GetJacobiDet();
              rhoi_shapes.Cols(2*i, 2*i+2) = shapes.Cols(2*i, 2*i+2) * rhoi;
            }
          elmat = rhoi_shapes * Trans(shapes);
          
          for (int i = 0; i < dofs.Size(); i++)
            for (int j = 0; j < dofs.Size(); j++)
              if (fabs(elmat(i,j)) > 1e-8)
                {
                  rowind.Append(dofs[i]);
                  colind.Append(dofs[j]);
                  values.Append(elmat(i,j));
                }
        }
      tint.Stop();

      tsp.Start();
      auto spmat = SparseMatrix<double>::CreateFromCOO(rowind, colind, values, GetNDof(), GetNDof());
      tsp.Stop();
      
      return make_shared<MySuperSparseMatrix> (std::move(*spmat));
    }
    else
      throw Exception("only 2d and 3d implemented");
    }


    
    std::map<ELEMENT_TYPE, IntegrationRule> 
      
      HDivPrimalCells::
      GetIntegrationRules(bool fix_lo) const
    {
      std::map<ELEMENT_TYPE, IntegrationRule> rules;
      rules[ET_TRIG] = PrimalCellIR(IR,false);

      rules[ET_TET] = PrimalVolIR(IR,false);

      //auto irseg = SelectIntegrationRule(ET_SEGM, 2*order+2);
      rules[ET_SEGM] = PrimalSegmIR(IR, false);

      if (IR.Size()==1 && fix_lo)
      {
        for (auto & ip : rules[ET_TRIG])
          ip.SetWeight(1.0/6);
        for (auto & ip : rules[ET_TET])
          ip.SetWeight(1.0/24);
      }


      return rules;
    }    
};    
