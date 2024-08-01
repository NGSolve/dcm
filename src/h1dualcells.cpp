#include <comp.hpp>
#include "h1dualcells.hpp"
#include "hdivprimalcells.hpp"
#include "intrules.hpp"
#include "supersparse.hpp"


namespace ngcomp
{

  // getting rid of the problems with ConstEBE<T> and make_shared
  //
  typedef decltype (ConstantElementByElementMatrix (5,5,Matrix<>(),
                                                  declval<Table<int>>(), declval<Table<int>>())) T_ConstEBE;

  class H1DualCellSegm : public ScalarFiniteElement<1>, public VertexOrientedFE<ET_SEGM>
  {
    const IntegrationRule & GaussRadauIR;
  public:
    H1DualCellSegm (int order, const IntegrationRule & _GaussRadauIR)
      : ScalarFiniteElement<1> (2*(order+1), order), GaussRadauIR(_GaussRadauIR)
    { ; }
    using VertexOrientedFE<ET_SEGM>::SetVertexNumbers;
    virtual ELEMENT_TYPE ElementType() const { return ET_SEGM; }

    
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceVector<> shape) const
    {
      ArrayMem<double, 20> polx(order+1);
      
      double lam[] = { ip(0), 1-ip(0) };
      int maxlam = 0;
      if (lam[1] > lam[0])
          maxlam = 1;

      shape.Range(GetNDof()) = 0;

      double lx = lam[(maxlam+1)%2];
      double xi = 2*lx;
      
      LagrangePolynomials(xi, GaussRadauIR, polx);
      shape(maxlam) = polx[0];
      
      int ii = 2;
      IVec<2> e = GetVertexOrientedEdge(0);
      for (int j = 0; j < 2; j++)
        {
          if (e[j] == maxlam)
            {
              for (int k = 0; k < order; k++)
                shape(ii+k) = polx[k+1];
            }
          ii+=order;
        }
    }

    virtual void CalcDShape (const IntegrationPoint & ip, 
                             BareSliceMatrix<> baredshape) const
    {
      ArrayMem<AutoDiff<1>, 20> polx(order+1);

      auto dshape = baredshape.AddSize(GetNDof(), 1);
      dshape = 0;
      
      double lam[] = { ip(0), 1-ip(0) };
      int maxlam = 0;
      if (lam[1] > lam[0])
          maxlam = 1;

      double lx = lam[(maxlam+1)%2];
      double xi = 2*lx;
      
      AutoDiff<1> adxi(xi, 0);
      LagrangePolynomials(adxi, GaussRadauIR, polx);
      double sig = 1;
      if (maxlam == 0)
        sig = -1;
      dshape(maxlam,0) = 2*sig*polx[0].DValue(0);
      
      int ii = 2;
      IVec<2> e = GetVertexOrientedEdge(0);
      for (int j = 0; j < 2; j++)
        {
          if (e[j] == maxlam)
            {
              for (int k = 0; k < order; k++)
                dshape(ii+k,0) = 2*sig*polx[k+1].DValue(0);
            }
          ii+=order;
        }
    }
  };




  
  class H1DualCellTrig : public ScalarFiniteElement<2>, public VertexOrientedFE<ET_TRIG>
  {
    const IntegrationRule & GaussRadauIR;
  public:
    H1DualCellTrig (int order, const IntegrationRule & _GaussRadauIR)
      : ScalarFiniteElement<2> (3*sqr(order+1), order), GaussRadauIR(_GaussRadauIR)
    { ; }
    using VertexOrientedFE<ET_TRIG>::SetVertexNumbers;
    virtual ELEMENT_TYPE ElementType() const { return ET_TRIG; }

    virtual void Interpolate (const ElementTransformation & trafo, 
                              const class CoefficientFunction & func, SliceMatrix<> coefs,
                              LocalHeap & lh) const
    {
      cout << "Interpolate called, ndof = " << ndof << endl;
      Matrix mat(ndof, ndof);
      Vector rhs(ndof);

      mat = 0.0;
      rhs = 0.0;
      
      auto [pnts, inds, nump] = GetNanoPoints<ET_TRIG> (order+1, true);
      cout << pnts.Size() << nump << endl;

      Vector shape(ndof);

      
      for (int i = 0; i < 3*(order+1)*(order+1); i++) // sum over real edges, only 
        {
          Vec<2> p = pnts[i];

          IntegrationPoint ip2d(p(0), p(1), 0,0);
          MappedIntegrationPoint<2,2> mip(ip2d, trafo);
          
          double fi = func.Evaluate (mip);

          CalcShape (ip2d, shape);

          mat.Row(i) = shape;
          rhs(i) = fi;
        }
      //cout << mat << endl;
      //cout << rhs << endl;
      CalcInverse (mat);
      coefs.Col(0) = mat * rhs;
    }
    
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceVector<> shape) const
    {
      ArrayMem<double, 20> polx(order+1), poly(order+1);      
      
      double lam[] = { ip(0), ip(1), 1-ip(0)-ip(1) };
      int maxlam = 0;
      for (int j = 1; j < 3; j++)
        if (lam[j] > lam[maxlam])
          maxlam = j;

      shape.Range(GetNDof()) = 0;

      double lx = lam[(maxlam+1)%3];
      double ly = lam[(maxlam+2)%3];
      double xi  = (3+2*lx-2*ly)/2 - sqrt( sqr( (2*ly-2*lx-3)/2 ) - 6*lx );      
      double eta = (3+2*ly-2*lx)/2 - sqrt( sqr( (2*lx-2*ly-3)/2 ) - 6*ly );
      
      LagrangePolynomials(xi, GaussRadauIR, polx);
      LagrangePolynomials(eta, GaussRadauIR, poly);
      shape(maxlam) = polx[0]*poly[0];
      
      int ii = 3;
      for (int i = 0; i < 3; i++)
        {
          IVec<2> e = GetVertexOrientedEdge(i);
          for (int j = 0; j < 2; j++)
            {
              if (e[j] == maxlam)
                {
                  if (e[1-j] == (maxlam+1)%3)
                    for (int k = 0; k < order; k++)
                      shape(ii+k) = polx[k+1]*poly[0];
                  else
                    for (int k = 0; k < order; k++)
                      shape(ii+k) = poly[k+1]*polx[0];
                }
              ii+=order;
            }
        }

      ii += maxlam*sqr(order);
      for (int i = 0; i < order; i++)
        for (int j = 0; j < order; j++)
          shape(ii++) = polx[i+1]*poly[j+1];
    }
    
    virtual void CalcDShape (const IntegrationPoint & ip, 
                             BareSliceMatrix<> baredshape) const
    {
      ArrayMem<AutoDiff<1>, 20> polx(order+1), poly(order+1);      
      
      double lam[] = { ip(0), ip(1), 1-ip(0)-ip(1) };
      int maxlam = 0;
      for (int j = 1; j < 3; j++)
        if (lam[j] > lam[maxlam])
          maxlam = j;

      auto dshape = baredshape.AddSize(GetNDof(), 2);
      dshape = 0;

      double lx = lam[(maxlam+1)%3];
      double ly = lam[(maxlam+2)%3];
      double xi  = (3+2*lx-2*ly)/2 - sqrt( sqr( (2*ly-2*lx-3)/2 ) - 6*lx );      
      double eta = (3+2*ly-2*lx)/2 - sqrt( sqr( (2*lx-2*ly-3)/2 ) - 6*ly );
      AutoDiff<1> adxi(xi, 0);
      AutoDiff<1> adeta(eta, 0);
      

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
      Mat<2,2> trafo = Trans(Inv(dxy_dxieta));
      
      LagrangePolynomials(adxi, GaussRadauIR, polx);
      LagrangePolynomials(adeta, GaussRadauIR, poly);

      auto assign_dshape = [trafo, dshape, &polx, &poly](int nr, int ix, int iy)
        {
          dshape.Row(nr) = trafo*Vec<2>(polx[ix].DValue(0)*poly[iy].Value(),
                                        polx[ix].Value()*poly[iy].DValue(0));
        };

      assign_dshape(maxlam, 0, 0);
      
      int ii = 3;
      for (int i = 0; i < 3; i++)
        {
          IVec<2> e = GetVertexOrientedEdge(i);
          for (int j = 0; j < 2; j++)
            {
              if (e[j] == maxlam)
                {
                  if (e[1-j] == (maxlam+1)%3)
                    for (int k = 0; k < order; k++)
                      assign_dshape(ii+k, k+1, 0);
                  else
                    for (int k = 0; k < order; k++)
                      assign_dshape(ii+k, 0, k+1);
                }
              ii+=order;
            }
        }

      ii += maxlam*sqr(order);
      for (int i = 0; i < order; i++)
        for (int j = 0; j < order; j++)
          assign_dshape(ii++, i+1, j+1);
    } 
  };



    H1DualCells::
      H1DualCells (shared_ptr<MeshAccess> ama, const Flags & flags)
      : FESpace (ama, flags)      
    {
    switch (ma->GetDimension())
    {
    case 1:
      {
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<1>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<1>>>();
        
        break;
      }
    case 2:
      {
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();
        break;
      }
    }


      GaussRadauIR = IntegrationRule();
      Array<double> xi, wi;
      ComputeGaussRadauRule (order+1, xi, wi);
      xi[0] += 1e-12;
      for (auto i : Range(xi))
        GaussRadauIR.Append (IntegrationPoint (xi[i], 0, 0, wi[i]));
      
      // cout << "GR-IR = " << GaussRadauIR << endl;
    }


    void H1DualCells::
      Update()
    {
      FESpace::Update();
      size_t ndof = ma->GetNV();
      
      first_edge_dofs.SetSize(ma->GetNEdges()+1);
      first_edge_dofs[0] = ndof;
      for (auto i : Range(ma->GetNEdges()))
        first_edge_dofs[i+1] = ndof += 2*order;
      if (ma->GetDimension()>1)
      {
        first_face_dofs.SetSize(ma->GetNFaces()+1);
        first_face_dofs[0] = ndof;
        for (auto i : Range(ma->GetNFaces()))
          first_face_dofs[i+1] = ndof += 3*sqr(order);
      }

      SetNDof(ndof);
    }
    
    void H1DualCells::
      GetDofNrs (ElementId ei, Array<DofId> & dnums) const
    {
      dnums.SetSize0();

      auto el = ma->GetElement(ei);
      for (auto v : el.Vertices())
        dnums.Append(v);
      for (auto ed : el.Edges())
        dnums += Range(first_edge_dofs[ed], first_edge_dofs[ed+1]);
      for (auto fa : el.Faces())
        dnums += Range(first_face_dofs[fa], first_face_dofs[fa+1]);
    }


    FiniteElement & H1DualCells::
      GetFE (ElementId ei, Allocator & alloc) const
    {
      auto ngel = ma->GetElement (ei);
      switch (ngel.GetType())
        {
        case ET_SEGM:
          {
            auto segm = new (alloc) H1DualCellSegm(order, GaussRadauIR);
            segm->SetVertexNumbers (ngel.vertices);
            return *segm;
          }

        case ET_TRIG:
          {
            auto trig = new (alloc) H1DualCellTrig(order, GaussRadauIR);
            trig->SetVertexNumbers (ngel.vertices);
            return *trig;
          }
        default:
          throw Exception("element not implemented");
        }
    }

  shared_ptr<BaseMatrix> H1DualCells::GetGradientOperator2D (bool dual) const
  {
    //throw Exception(" does not work yet...");
    cout << "called GetGradientOperator" << endl;
    Flags flags;
    flags.SetFlag ("order", order-1);
    flags.SetFlag ("dirichlet", this->flags.GetNumListFlag("dirichlet"));
    flags.SetFlag ("uniform_order", true);
    auto feshdiv = make_shared<HDivPrimalCells>(ma, flags);
    feshdiv -> Update();
    feshdiv -> FinalizeUpdate();
    LocalHeap lh(10*1000*1000);
    shared_ptr<BaseMatrix> sum;
    
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
    
    IntegrationRule irsegm(ET_SEGM, 2*order+10);

    auto altshapes = true;
    
    auto [pnts,ind,num] = GetNanoPoints<ET_TRIG> (order+1);
    auto [lines,line_ends] = GetNanoLines<ET_TRIG> (order+1);

    for (auto elclass_inds : table)
      {
        if (elclass_inds.Size() == 0) continue;
        
        ElementId ei(VOL,elclass_inds[0]);
        
        auto & felh1 = dynamic_cast<const ScalarFiniteElement<2>&> (GetFE (ei, lh));
        auto & felhdiv = dynamic_cast<const HDivPrimalCellTrig&> (feshdiv->GetFE (ei, lh));
        auto & trafo = ma->GetTrafo(ei, lh);

        //cout << felhdiv.GetNDof() << endl; 
        Matrix<> gradmat(felhdiv.GetNDof(), felh1.GetNDof());
        Matrix<> mixedmass(felhdiv.GetNDof(), felhdiv.GetNDof());
        Matrix<> lumpedmass(felhdiv.GetNDof(), felhdiv.GetNDof());
        
        gradmat = 0.;
        mixedmass = 0.;
        lumpedmass = 0.;

        Matrix shapehdiv(felhdiv.GetNDof(), 2);        
        Vector shapeh1(felh1.GetNDof());
        
        bool altshapeshdiv = true;


          // gradmat
          for (auto i : Range(line_ends))
          {
            IntegrationPoint p(pnts[line_ends[i][0]][0], pnts[line_ends[i][0]][1],0);
            felh1.CalcShape(p, shapeh1);    
            gradmat.Row(i) += -shapeh1;
            IntegrationPoint q(pnts[line_ends[i][1]][0], pnts[line_ends[i][1]][1],0);
            felh1.CalcShape(q, shapeh1);   
            gradmat.Row(i) += shapeh1;
          }
        //cout << "gradmat = " << gradmat << endl;
        //throw Exception("test");
      
          // mixedmass
          for (auto i : Range(lines))
          {
            //cout << i << endl;
            for (auto e : lines[i])
              {
                //cout << e[0] << endl;
                Vec<2> micropts[2] =
                  { pnts[e[0]], pnts[e[1]]};
                //cout << micropts[1] << micropts[0]<< endl;
                Vec<2> tau = micropts[1]-micropts[0];
            
                  for (auto ip : irsegm)
                    {
                      Vec<2> p = (1-ip(0)) * micropts[0] + ip(0) * micropts[1];
                      IntegrationPoint ipseg(p(0), p(1), 0);
                      if (altshapeshdiv)
                        felhdiv.CalcAltShape(ipseg, shapehdiv);    
                      else
                        felhdiv.CalcShape(ipseg, shapehdiv);    
                    mixedmass.Row(i) += ip.Weight()*(shapehdiv*tau);
                    }
                  }
          }
              
        //cout << "mixedmass = " << mixedmass << endl;
        CalcInverse (mixedmass);
        
        //cout << "mixedmassinv = " << mixedmass << endl;
            // integrate
            // lumpedmass ...
            //lumped integration is the same (except order=0 fixup)
            //auto irs = feshdiv->GetIntegrationRules();
            //auto irs = GetIntegrationRules(2*order+6);
            //IntegrationRule ir = move(irs[ET_TRIG]);
            
            auto irs = ngcomp::GetIntegrationRules(2*order+6);
            //auto irs = feshdiv->GetIntegrationRules();
            IntegrationRule ir = move(irs[ET_TRIG]);            
            
            Matrix shapes(felhdiv.GetNDof(), ir.Size()*2);
            Matrix Piolashapes(felhdiv.GetNDof(), ir.Size()*2);
        
            for (int i = 0; i < ir.Size(); i++)
              {
                felhdiv.CalcAltShape (ir[i], shapes.Cols(2*i, 2*i+2));
                felhdiv.CalcAltShape (ir[i], Piolashapes.Cols(2*i, 2*i+2));
                //felhdiv.CalcShape (ir[i], shapes.Cols(3*i, 3*i+3));
                //felhdiv.CalcPiolaShape (ir[i], Piolashapes.Cols(3*i, 3*i+3));
                Piolashapes.Cols(2*i, 2*i+2) *= ir[i].Weight();
              }

            lumpedmass = shapes * Trans (Piolashapes);
        
            //cout << "lumpedmass = " << lumpedmass << endl;

        Matrix<> tmpmat(felhdiv.GetNDof(), felh1.GetNDof());
        tmpmat = mixedmass * gradmat;

        if (dual)
          gradmat = lumpedmass * tmpmat;
        else
          gradmat = tmpmat;
        //cout << "finalmat = " << finalmat << endl;
        
        Table<DofId> xdofs(elclass_inds.Size(), felh1.GetNDof()),
          ydofs(elclass_inds.Size(), felhdiv.GetNDof());
        
        Array<DofId> dnumsx, dnumsy;
        for (auto i : Range(elclass_inds))
          {
            ElementId ei(VOL, elclass_inds[i]);
            feshdiv->GetDofNrs(ei, dnumsy);
            GetDofNrs(ei, dnumsx);
            xdofs[i] = dnumsx;
            ydofs[i] = dnumsy;
          }
        auto mat = make_shared<T_ConstEBE>
          ( feshdiv->GetNDof(),GetNDof(),
           gradmat, std::move(ydofs),std::move(xdofs));
        if (sum)
          sum = make_shared<SumMatrix>(sum, mat);
        else
          sum = mat;
      }
    return sum;
    
  }
  
  shared_ptr<BaseMatrix> H1DualCells ::
  GetMassOperator(shared_ptr<CoefficientFunction> rho, shared_ptr<Region> defon, LocalHeap & lh) const
  {
    HeapReset hr(lh);
    auto irs = GetIntegrationRules(nullopt);
    IntegrationRule ir = move(irs[ET_TRIG]);

    auto & felref = dynamic_cast<const H1DualCellTrig&> (GetFE(ElementId(VOL,0), lh));
    Matrix shapes(felref.GetNDof(), ir.Size());
    Matrix rhoi_shapes(felref.GetNDof(), ir.Size());
    
    for (int i = 0; i < ir.Size(); i++)
      felref.CalcShape (ir[i], shapes.Col(i));

    for (int i = 0; i < shapes.Height(); i++)
      for (int j = 0; j < shapes.Width(); j++)
        if (fabs(shapes(i,j)) < 1e-8)
          shapes(i,j) = 0;
    rhoi_shapes = shapes;
    
    Array<int> rowind, colind;
    Array<double> values;
    Matrix<> elmat(felref.GetNDof());
    Array<DofId> dofs(felref.GetNDof());
    
    for (size_t nr : Range(ma->GetNE()))
      {
        if (defon && !defon->Mask()[ma->GetElIndex(ElementId(VOL,nr))]) continue;
        
        HeapReset hr(lh);
        MappedIntegrationRule<2,2> mir(ir, ma->GetTrafo(ElementId(VOL,nr), lh), lh);
        GetDofNrs(ElementId(VOL,nr), dofs);

        FlatVector<> rhovals(mir.Size(), lh);
        rhovals = 1;
        rho->Evaluate (mir, rhovals.AsMatrix(mir.Size(), 1));
        
        for (size_t i = 0; i < mir.Size(); i++)
          {
            double rhoi = rhovals(i);
              
            //Mat<2,2> F = mir[i].GetJacobian();
            rhoi *= ir[i].Weight() * mir[i].GetJacobiDet();
            rhoi_shapes.Col(i) = rhoi * shapes.Col(i);
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

    auto spmat = SparseMatrix<double>::CreateFromCOO(rowind, colind, values, GetNDof(), GetNDof());
    
    return make_shared<MySuperSparseMatrix> (move(*spmat));
  }



  
    std::map<ELEMENT_TYPE, IntegrationRule> H1DualCells::
      GetIntegrationRules(optional<int> intorder) const
    {
      std::map<ELEMENT_TYPE, IntegrationRule> rules;

      IntegrationRule IR(GaussRadauIR);
      if (intorder.has_value())
      {
        IR = IntegrationRule();
        Array<double> xi, wi;
        ComputeGaussRadauRule (intorder.value(), xi, wi);
        xi[0] += 1e-12;
        for (auto i : Range(xi))
          IR.Append (IntegrationPoint (xi[i], 0, 0, wi[i]));
      }
      rules[ET_TRIG] = PrimalCellIR(IR);

      // auto irseg = SelectIntegrationRule(ET_SEGM, 2*order);
      //auto irseg = SelectIntegrationRule(ET_SEGM, intorder.value_or(2*order+2));
      
      //rules[ET_SEGM] = PrimalSegmIR(irseg);

      IntegrationRule irseg;
      for (auto & ip : GaussRadauIR)
        {
          double x = ip(0), w = ip.Weight();
          irseg.Append ( IntegrationPoint(x/2+1e-12, 0, 0, 0.5*w));
          irseg.Append ( IntegrationPoint(1-x/2-1e-12, 0, 0, 0.5*w));
        }
      rules[ET_SEGM] = move(irseg);
      
      return rules;
    }    
    


}
