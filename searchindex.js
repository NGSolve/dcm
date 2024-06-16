Search.setIndex({"docnames": ["examples", "examples/qualitative_dispersion", "examples/ring_resonator", "examples/waveguide_3d", "installation", "intro", "maxwell", "maxwell/diffops", "maxwell/hcurl", "maxwell/introduction", "maxwell/mass_lumping", "maxwell/time_stepping"], "filenames": ["examples.md", "examples/qualitative_dispersion.ipynb", "examples/ring_resonator.ipynb", "examples/waveguide_3d.ipynb", "installation.md", "intro.md", "maxwell.md", "maxwell/diffops.md", "maxwell/hcurl.md", "maxwell/introduction.md", "maxwell/mass_lumping.md", "maxwell/time_stepping.md"], "titles": ["<span class=\"section-number\">3. </span>Examples", "<span class=\"section-number\">3.1. </span>A plane wave on a square", "<span class=\"section-number\">3.2. </span>Ring resonator", "<span class=\"section-number\">3.3. </span>Dualmesh method for an electromagnetic waveguide in 3d", "<span class=\"section-number\">1. </span>Installation", "The Dual Cell Method in NGSolve", "<span class=\"section-number\">2. </span>The dual cell method for the time-domain Maxwell system", "<span class=\"section-number\">2.5. </span>The discrete differential operators", "<span class=\"section-number\">2.3. </span>Discrete spaces on dual cells", "<span class=\"section-number\">2.2. </span>Introduction", "<span class=\"section-number\">2.4. </span>Mass lumping", "<span class=\"section-number\">2.6. </span>Time stepping"], "terms": {"A": [0, 5, 8, 9], "plane": [0, 3, 5], "wave": [0, 2, 5, 7, 9], "squar": [0, 5, 8, 10], "ring": [0, 5], "reson": [0, 5], "we": [1, 2, 4, 5, 6, 7, 8, 9, 10, 11], "solv": [1, 2, 7, 9], "two": [1, 2, 5, 8, 9], "dimension": [1, 8, 9], "equat": [1, 2, 5, 6, 11], "find": [1, 6, 7, 8, 9], "h": [1, 2, 3, 6, 7, 8, 9], "0": [1, 2, 3, 6, 7, 8, 9, 10, 11], "t": [1, 2, 3, 6, 7, 9], "1": [1, 2, 3, 7, 8, 9, 10, 11], "omega": [1, 2, 3, 6, 8, 9], "vector": [1, 2, 7, 8], "field": [1, 2, 5, 6, 9], "e": [1, 2, 3, 4, 5, 6, 7, 8, 9], "mathrm": [1, 2, 6, 7, 9], "div": [1, 2, 7], "ar": [1, 2, 5, 6, 7, 8, 9, 10, 11], "begin": [1, 2, 6, 11], "align": [1, 2, 6, 11], "partial_t": [1, 2, 6, 8, 9, 11], "x": [1, 3, 6, 7, 8], "nabla": [1, 2], "f": [1, 2, 7, 8], "exp": [1, 2, 3, 7], "400": [1, 2], "y": [1, 2, 3, 7], "2": [1, 2, 3, 5, 7, 8, 9, 10, 11], "partial": [1, 7, 8, 9], "end": [1, 2, 3, 6, 8, 11], "from": [1, 2, 3, 4, 7, 8, 10], "ngsolv": [1, 2, 3, 7, 8, 10], "import": [1, 2, 3, 7, 8, 10], "dualcellspac": [1, 2, 3, 7, 8, 10], "dc": [1, 2, 3, 7, 8, 10], "time": [1, 3, 5, 7, 8, 9, 10], "webgui": [1, 2, 3, 7, 8], "draw": [1, 2, 3, 7, 8], "after": 1, "necessari": [1, 9], "defin": [1, 2, 7, 8, 9, 10], "some": [1, 6, 7], "paramet": [1, 2, 3, 9], "mesh": [1, 2, 3, 5, 7, 9, 10], "maxh": [1, 2, 3, 7, 8, 10], "03": [1, 2, 7], "tend": [1, 2, 3], "order": [1, 2, 3, 4, 5, 7, 8, 10], "h0": 1, "cf": [1, 2, 3, 7], "20": [1, 2, 5], "e0": 1, "unit_squar": [1, 7, 8, 10], "generatemesh": [1, 2, 3, 7, 8, 10], "space": [1, 7, 10], "fesh": [1, 7], "h1dualcel": [1, 2, 7], "fese": [1, 7], "hdivprimalcel": [1, 2, 7], "To": [1, 2, 6, 8, 10], "mass": [1, 2, 3, 5, 7, 9, 11], "bilinear": [1, 3, 8], "form": [1, 3, 9, 11], "need": [1, 2, 6, 7, 10], "special": [1, 3], "integr": [1, 2, 3, 7, 8, 9, 10], "rule": [1, 2, 3, 7, 10], "de": [1, 3, 7], "tnt": [1, 2, 3, 7, 10], "dh": [1, 3], "dxh": [1, 3], "dx": [1, 2, 3, 7, 10], "intrul": [1, 2, 3, 7, 10], "getintegrationrul": [1, 2, 3, 7, 10], "dsw": [1, 2, 3, 7], "element_boundari": [1, 2, 3, 7], "true": [1, 2, 3, 7, 8], "6": [1, 2, 8], "dxw": [1, 2, 3, 7], "id": [1, 2], "massh": 1, "massinv": 1, "invers": [1, 2, 3, 7, 9, 10], "massinvh": 1, "normal": [1, 2, 3, 7], "specialcf": [1, 2, 3, 7], "grad": [1, 7], "bilinearform": [1, 2, 3, 7, 10], "geom_fre": [1, 2, 3, 7], "assembl": [1, 2, 3, 7, 10], "mat": [1, 2, 3, 7, 10], "lffh": 1, "linearform": [1, 2, 3, 7], "The": [1, 2, 4, 9, 10], "maxim": 1, "admiss": [1, 9], "step": [1, 3, 8, 9, 10], "mai": [1, 2, 4, 6, 8, 9, 10], "estim": 1, "us": [1, 2, 3, 7, 8, 9, 10, 11], "simpl": [1, 4, 8], "power": 1, "iter": 1, "def": [1, 2], "estimate_tau": 1, "maxstep": 1, "1000": 1, "tol": 1, "1e": [1, 3, 10], "4": [1, 2, 3, 7, 8, 11], "vec": [1, 2, 3, 7, 8], "createcolvector": 1, "setrandom": [1, 7, 10], "tmp": [1, 2, 10], "createvector": [1, 2, 3, 10], "lam": 1, "i": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], "rang": [1, 7, 8, 10], "print": [1, 2, 3, 7, 8, 10], "r": [1, 2, 3, 6, 7, 8, 9], "data": [1, 2, 3, 7, 10], "lamnew": 1, "innerproduct": [1, 3], "tau": [1, 2, 3, 9, 11], "sqrt": [1, 2], "re": 1, "norm": 1, "diff": 1, "return": [1, 2], "did": 1, "converg": 1, "last": 1, "timestep": [1, 2, 3, 9, 11], "format": [1, 2, 3, 7, 8, 10], "9": [1, 2], "438193e": 1, "It": [1, 2, 5, 9, 11], "remain": [1, 2, 9], "set": [1, 2, 3, 8], "initi": 1, "condit": [1, 6], "draweveri": [1, 2, 3], "30": 2, "gfe": [1, 3, 7], "gridfunct": [1, 2, 3, 7, 8], "gfh": [1, 3, 7], "gfh_histori": 1, "multidim": [1, 2, 3, 8], "scene": [1, 2, 3], "intpoint": [1, 2, 8], "getwebguipoint": [1, 2, 3, 7, 8], "autoscal": [1, 2, 3], "fals": [1, 2, 3], "min": [1, 2, 3, 8], "max": [1, 2, 3, 8], "start": [1, 2, 8, 9, 10], "loop": 1, "now": [1, 2, 3, 7, 10], "nowstart": 1, "energi": [1, 9], "tmph": 1, "tmpe": [1, 3], "subtim": 1, "taskmanag": [1, 2, 3, 7, 10], "while": [1, 2, 3, 9], "timepass": [1, 3], "before_energy_tim": 1, "addmultidimcompon": [1, 2, 3], "redraw": [1, 2, 3], "append": [1, 2], "current": 1, "dof": [1, 2, 3, 7, 8], "": [1, 2, 3, 7, 8, 9, 10], "ndof": [1, 2, 3, 7, 8, 10], "comptim": 1, "n": [1, 2, 6, 7, 8, 10], "per": [1, 8], "second": [1, 2], "1545": 1, "1491700646156": 1, "162256": 1, "653121e": [], "07": [2, 3], "anim": [1, 2, 3, 8], "observ": [1, 10], "preserv": 1, "modifi": [1, 4], "discret": [1, 2, 9, 11], "matplotlib": [1, 2, 3, 4, 10], "pyplot": [1, 2, 3, 10], "pl": [1, 2, 3, 10], "plot": [1, 2, 3, 10], "ylim": 1, "simul": [2, 5], "propag": 2, "an": [0, 2, 4, 5, 8], "electromagnet": [0, 2, 5], "tm": 2, "through": 2, "devic": 2, "dual": [2, 10], "cell": [1, 2, 7, 10], "method": [0, 1, 2, 10], "add": 2, "desir": 2, "sketch": 2, "below": 2, "where": [2, 5, 6, 8, 9, 10, 11], "horizont": 2, "waveguid": [0, 2, 5], "suppos": 2, "unbound": 2, "govern": 2, "scalar": [2, 9], "p": [2, 7, 8, 9, 10, 11], "u": [2, 7, 10], "l": [2, 5, 8, 9], "int_": [2, 6, 7], "cdot": [2, 7, 8, 9], "v": [2, 7, 9, 10, 11], "pq": 2, "q": [2, 7, 8, 9], "fq": 2, "all": [2, 6, 7, 8], "test": 2, "function": [1, 2, 5, 7, 9, 10], "domain": [2, 5, 7, 8, 9], "surround": 2, "absorb": 2, "layer": 2, "damp": 2, "system": [2, 5, 9, 11], "alpha": 2, "omega_": 2, "tb": 2, "cup": 2, "lr": 2, "top": [2, 7, 8, 11], "mathbf": [2, 6, 7, 8, 9, 11], "hat": [2, 8, 9], "left": [2, 8, 9, 11], "nn": 2, "right": [2, 7, 8, 9, 11], "int": [2, 8], "fv": 2, "omega_c": 2, "direct": [2, 8], "case": 2, "given": [2, 8, 9, 11], "parallel": 2, "wire": 2, "shape": [2, 10], "one": [2, 6, 7, 8, 9, 10], "between": [2, 9, 10], "splinegeometri": 2, "netgen": [2, 3], "geom2d": 2, "gm": 2, "numpi": [2, 3, 4, 10], "np": [2, 3, 10], "geo": [2, 3], "xneg": 2, "43": 2, "xpo": 2, "yneg": 2, "48": [2, 8], "ypo": 2, "wslab": 2, "04": 2, "cringx": 2, "cringi": 2, "rring": 2, "gap": 2, "005": 2, "pntx": 2, "pnty": 2, "pt": 2, "yi": 2, "xi": [2, 8], "addpoint": 2, "inner": [2, 8, 9], "rect": 2, "line": [2, 9], "leftdomain": 2, "rightdomain": 2, "3": [2, 3, 6, 7, 8, 9, 10], "5": [2, 3, 7, 8, 10], "bc": 2, "normal_wg_rightbottom": 2, "normal_wg_leftbottom": 2, "7": [2, 8], "normal_wg_righttop": 2, "8": [2, 3], "normal_wg_lefttop": 2, "11": 2, "10": [2, 3, 5, 10], "addcircl": 2, "c": [2, 3, 6, 8, 9], "setmateri": 2, "air": 2, "eps_nin": 2, "thi": [2, 3, 5, 6, 7, 8, 9, 10, 11], "result": [2, 7, 10], "follow": [2, 7, 8, 9], "triangul": 2, "mesh_inn": 2, "curv": [2, 3], "howev": [2, 9], "also": [2, 4, 6, 9], "want": [2, 4], "perfectli": 2, "match": 2, "pmlwidth": 2, "05": 2, "createpml": 2, "ha": [2, 4, 8, 9], "ad": 2, "our": [2, 7, 8, 9, 10], "interior": 2, "bd": 2, "getboundari": 2, "default": 2, "getmateri": 2, "pml_default": 2, "pml_corner": 2, "pml_default_duplicate_1": 2, "pml_normal_wg_rightbottom": 2, "pml_default_duplicate_2": 2, "pml_normal_wg_righttop": 2, "pml_default_duplicate_3": 2, "pml_default_duplicate_4": 2, "pml_default_duplicate_5": 2, "pml_normal_wg_lefttop": 2, "pml_default_duplicate_6": 2, "pml_normal_wg_leftbottom": 2, "pml_default_duplicate_7": 2, "boundari": [1, 2, 3, 5, 6, 7, 8, 9], "wavelength": 2, "542": 2, "fcen": 2, "tpeak": 2, "sourcei": 2, "coefficientfunct": 2, "sin": [2, 3], "pi": [2, 3], "fes_facet": 2, "facetfespac": 2, "gfsourc": 2, "definedon": [2, 3], "t_envelop": 2, "ab": 2, "els": 2, "delta": 2, "001": 2, "arang": [2, 3], "xlabel": 2, "figur": [2, 10], "materi": 2, "background": 2, "medium": [2, 6], "eps_r": 2, "startswith": 2, "pml_normal_wg": 2, "ep": 2, "scale": [2, 9], "non": 2, "corner": 2, "nvec": 2, "cfn": 2, "next": [2, 8], "up": [2, 7, 8], "fes_u": 2, "fes_p": 2, "dirichlet": [2, 3], "outer": 2, "vectori": [2, 8, 9, 10], "fe": [2, 3, 10], "total": [2, 3], "030700e": [], "887000e": [], "383540e": [], "which": [2, 7, 8, 9, 10], "respect": [2, 5, 8, 9, 10], "These": 2, "can": [2, 4, 8, 10], "obtain": [2, 7, 8, 9, 10], "via": [2, 4, 8, 10], "over": [1, 2, 3, 7], "volum": 2, "element": [1, 2, 4, 5, 7, 8, 9, 10], "ir": [2, 7, 10], "et": [2, 10], "segm": 2, "fem": 2, "integrationrul": [2, 10], "object": [2, 3], "0x70232cb071f0": [], "trig": [2, 10], "0x70232cb06930": [], "tet": 2, "0x70232cb07630": [], "diverg": 2, "contain": [2, 8], "jump": 2, "term": [1, 2, 7, 9], "due": [2, 7, 9], "fact": [2, 9, 10], "v_h": 2, "discontinu": [2, 5, 9, 10], "b": [2, 5, 6, 7, 11], "matric": [2, 5, 7, 10, 11], "block": [1, 2, 5, 10], "diagon": [1, 2, 5, 10], "lump": [2, 5, 7], "thei": [2, 10], "implement": [2, 5, 7, 8, 9, 10], "effici": 2, "fespac": [2, 10], "invmassp": 2, "freedof": [2, 3], "invmassu": 2, "dim": 2, "pml1d": 2, "pml_normal": 2, "dampp1": 2, "dampp2": 2, "dampu1": 2, "outerproduct": 2, "dampu2": 2, "big": 2, "oper": [2, 3, 9, 11], "embed": 2, "small": 2, "emb_p": 2, "emb_phat": 2, "emb_u": 2, "emb_uhat": 2, "b_big": 2, "dampu_big": 2, "dampp_big": 2, "invmassp_big": 2, "invmassu_big": 2, "lastli": [2, 7], "q_big": [], "testfunct": 7, "lsrc": 2, "gfu": [2, 7, 8], "gfp_histori": 2, "compon": [2, 3, 9], "4e": [], "5e": [], "25": 8, "startnow": 2, "200": [], "mdofss": [], "087120000001544": [], "113300e": [], "ssss": [], "keyboardinterrupt": [], "traceback": [], "most": [4, 9], "recent": [4, 9], "call": 8, "In": [6, 8, 10], "37": [], "23": [], "21": [], "22": [], "24": [], "finish": [], "packag": [4, 8, 10], "high": [4, 5, 10], "finit": [4, 5, 8, 9, 10], "librari": [3, 4, 5], "thu": [4, 5, 7, 8, 9, 10], "main": [4, 9, 10], "premis": 4, "suffici": 4, "version": 4, "beforehand": 4, "wai": [4, 7, 8], "python": 4, "m": [4, 5, 10], "scipi": 4, "jupyt": 4, "ipyparallel": 4, "scikit": 4, "build": 4, "upgrad": 4, "webgui_jupyter_widget": 4, "For": [4, 5, 8, 9, 10, 11], "troubleshoot": 4, "refer": [4, 7, 8, 10], "variou": 4, "tutori": 4, "ngs24": 4, "document": 4, "If": 4, "you": 4, "do": 4, "updat": 4, "your": 4, "might": 4, "consid": 4, "virtual": 4, "environ": 4, "As": [4, 7, 10], "avail": 4, "core": 4, "again": [4, 8, 9], "have": [4, 8, 9, 10], "binari": 4, "git": 4, "http": [4, 5], "github": 4, "com": 4, "dcm": [4, 5], "built": 4, "sourc": [1, 4, 6], "dev": 4, "wa": [4, 10], "g": [4, 5, 9], "pre": 4, "pybind11_stubgen": 4, "isol": 4, "modul": 4, "clone": 4, "them": [4, 8], "done": [4, 7, 8, 9], "either": 4, "cmake": 4, "cd": 4, "mkdir": 4, "make": [4, 9], "j4": 4, "whether": [4, 10], "demo": 4, "dc_intrul": 4, "test_spac": 4, "wess": 5, "j": [5, 8, 9, 11], "sch\u00f6berl": 5, "tu": 5, "wien": 5, "institut": 5, "analysi": 5, "scientif": 5, "comput": [2, 5, 10], "base": 5, "joint": 5, "work": 5, "kapidani": 5, "codecasa": 5, "book": 5, "design": 5, "provid": [5, 10], "introduct": 5, "exampl": 5, "galerkin": 5, "acoustiv": 5, "mix": [3, 5], "formul": [5, 6, 8, 9], "disconitinu": 5, "variant": 5, "approxim": [5, 9, 10], "conform": [5, 8], "each": [5, 7, 8, 9, 10], "other": [5, 9], "ansatz": 5, "featur": [5, 9], "differ": [5, 8, 9], "full": [5, 7], "mathemat": 5, "kcschoberl21": 5, "wkcs23": [], "instal": 5, "maxwel": 5, "bernard": 5, "lorenzo": 5, "joachim": 5, "sch": 5, "\u00f6": 5, "berl": 5, "arbitrari": 5, "depend": [5, 8, 9, 10], "d": [3, 5, 6, 7], "axwel": 5, "phy": 5, "433": 5, "paper": 5, "No": 5, "110184": 5, "2021": 5, "doi": 5, "1016": 5, "jcp": 5, "marku": 5, "polynomi": [5, 8, 10], "degre": [5, 9], "acoust": 5, "2023": [], "arxiv": [], "2312": [], "14716": [], "section": [6, 7, 8, 9, 10], "describ": 6, "detail": 6, "how": 6, "appli": [6, 9, 10], "neglect": 6, "possibl": [6, 7, 9, 10], "here": [6, 10], "curl": [3, 6, 7], "text": [6, 9], "suitabl": [1, 6], "subset": 6, "mathbb": [6, 7, 8, 9], "close": 6, "constitut": 6, "relat": 6, "color": 6, "emph1": 6, "varepsilon": [6, 8, 9], "emph2": 6, "mu": [6, 8, 9], "permitt": 6, "permeabl": 6, "question": [6, 10], "weak": [6, 7, 8, 9], "eh": 6, "assum": 6, "homogen": 6, "abov": [6, 7, 8], "rewritten": 6, "377423e": [], "08": [1, 3], "270705e": [], "187766e": [], "version_info": [], "modulenotfounderror": [], "input": [], "name": 3, "9933350511855035": [], "1540": [], "06247129244141534": [], "381633e": [], "088": [], "160984e": [], "392659e": [], "35": [], "409019e": [], "322534e": [], "371604e": [], "0x7abf0dc794b0": [], "0x7abf0dc793f0": [], "0x7abf0dc78df0": [], "036663e": [], "3e": 2, "300": [], "dofss": 2, "when": [9, 10], "type": [], "problem": [7, 8, 10], "approach": [8, 9, 10], "typic": 9, "choic": [9, 10], "implicit": 9, "explicit": 10, "former": 9, "matrix": [7, 8, 9, 10, 11], "compos": [8, 9], "differenti": [1, 9, 11], "although": 9, "scheme": [9, 11], "uncondition": 9, "stabl": [3, 9, 11], "independ": 9, "factor": [9, 10], "larg": 9, "doe": 9, "well": [9, 11], "veri": [9, 10], "number": [8, 9], "freedom": 9, "hand": [7, 9], "mere": [7, 9], "There": 9, "exist": [8, 9], "sever": 9, "conveni": 9, "downsid": 9, "stabil": 9, "qualiti": 9, "finer": 9, "largest": 9, "guarante": 9, "popular": 9, "afdsj": [], "adf": [], "dg": 9, "appelogr20": [], "daniel": [], "appel": [], "fortino": [], "garcia": [], "olof": [], "runborg": [], "waveholtz": [], "solut": 3, "helmholtz": 9, "siam": [], "journal": 5, "42": [], "a1950": [], "a1983": [], "2020": 9, "1137": [], "19m1299062": [], "apart": [9, 10], "obviou": [9, 10], "better": 9, "applic": [7, 9, 10], "solver": 9, "interest": 9, "frequenc": 9, "scatter": [], "precondition": 9, "stolk": [], "eigenvalu": [9, 10], "lite": 9, "pinvit": 9, "lobpcg": 9, "earli": 9, "2000": 9, "knyazev": [], "filter": 9, "2e5": 9, "hit": 9, "googl": 9, "scholar": 9, "2e4": 9, "sinc": [8, 9, 10], "2022": 9, "low": 9, "onli": [8, 9], "hexahedr": [8, 9], "grid": [8, 9], "5e4": 9, "6e3": 9, "numer": [9, 10], "flux": 9, "penalti": 9, "both": 9, "been": 9, "around": 9, "late": 9, "1960": 9, "1970": 9, "still": 9, "wide": 9, "engin": 9, "sto21": [], "christiaan": [], "a3469": [], "a3502": [], "20m1359997": [], "kny01": [], "andrew": [], "toward": [], "optim": [1, 10], "precondit": [], "eigensolv": [], "local": [8, 9], "conjug": [], "gradient": [], "517": [], "541": [], "2001": [], "s1064827500366124": [], "nw24": [], "lothar": [], "nannen": [], "krylov": [], "2024": 5, "2402": [], "08515": [], "idea": [8, 10], "stoke": 9, "theorem": 9, "my": [], "caption": [], "visual": [3, 8, 9], "lead": 9, "quantiti": 9, "interlac": 9, "pointwis": 9, "mimmick": 9, "l_": [8, 9], "k": [8, 9], "point": [1, 3, 7, 8, 9, 10], "valu": 9, "x_i": 9, "y_j": 9, "z_k": 9, "tangenti": [8, 9], "stagger": 9, "h_": 9, "tfrac": [7, 9], "approx": 9, "e_": 9, "satisfi": 9, "leap": [9, 11], "frog": [9, 11], "size": 9, "_": [8, 9], "frac": [7, 9, 11], "quad": [8, 9, 11], "carri": 9, "goal": 9, "gener": 9, "tetrahedr": [8, 9], "157738e": [], "basi": [1, 7, 9, 10], "focuss": 9, "go": 9, "dimens": [9, 10], "primal": [1, 8, 9, 10], "interpret": 9, "orang": 9, "circl": 9, "piecewis": 9, "constant": 9, "dark": 9, "grei": 9, "consist": [8, 9], "four": [8, 9], "proce": 9, "unknown": 9, "neighbour": [8, 9], "continu": [8, 9], "tild": [7, 9], "c_j": 9, "bigcup_": [8, 9], "c_": 9, "mathcal": [7, 9], "_j": 9, "w_": 9, "setr": [], "foral": [8, 9], "cont": 9, "v_": [9, 11], "global": [7, 8, 9], "w": 9, "zero": [8, 9, 10], "outsid": [8, 9], "origin": 9, "pose": 9, "sum": [], "rot": 9, "sum_": [7, 8, 9], "semi": [9, 11], "ultra": 9, "variat": 9, "remark": 9, "aris": 9, "across": 9, "part": 9, "help": 9, "skew": 9, "symmetri": 9, "immedi": 9, "conserv": 9, "simplex": 9, "motiv": 9, "construct": 10, "symbol": [1, 3], "correspond": [1, 10], "nodal": [1, 8, 10], "higher": 1, "sparsiti": [1, 10], "pattern": [1, 10], "328608e": 1, "596200e": 2, "352700e": 2, "824640e": 2, "0x7efeaab2eef0": 2, "0x7efeaab2f5f0": 2, "0x7efeaab2c9f0": 2, "879820000000851": 2, "030279e": 2, "sss": 2, "averag": [2, 3], "1546": [], "424357e": [], "home": [], "mwess": [], "lib": [], "python3": [], "site": [], "project": 7, "__init__": [], "py": 10, "63": [], "userwarn": [], "unabl": [], "axes3d": [], "multipl": [], "being": [], "pip": [], "3d": [0, 5], "warn": 3, "0x737388a1c9b0": [], "0x737388a1c270": [], "0x73738c7f48b0": [], "11327": [], "613625e": [], "309617e": [], "0x7c895f45ba70": [], "0x7c895d73d2b0": [], "0x7c895f343330": [], "552361e": [], "50": 1, "342434e": [], "498972e": [], "265443e": [], "0x7f00a4f15030": [], "0x7f00a4fa5530": [], "0x7f00a4fa5d30": [], "735248e": [], "274404e": [], "0x71581d5d7330": [], "0x71581d5d4330": [], "0x71581d5d5df0": [], "586595e": [], "408633e": [], "0x7aa7f4fbe670": [], "0x7aa7f4fbe0b0": [], "0x7aa7f4fbff70": [], "drawtim": 2, "drawnow": 2, "computetim": 2, "45": 8, "34500765800476": [], "553796e": [], "wkcs24": 5, "physic": [5, 7, 8], "page": 5, "113196": 5, "org": 5, "decomposit": 8, "subdomain": 8, "bar": 8, "usual": [8, 10], "triangl": 8, "tetrahedra": 8, "therebi": 8, "barycent": 8, "connect": 8, "midpoint": 8, "edg": 8, "creat": 8, "quadrilater": 8, "share": 8, "vertex": 8, "triangular": 8, "split": 8, "becom": 8, "center": 8, "red": 8, "blue": 8, "uniqu": 8, "varphi_c": 8, "node": [8, 10], "x_0": [8, 10], "ldot": [8, 10], "unit": [8, 10], "interv": [8, 10], "tensor": [8, 10], "x_": 8, "xj": 8, "_p": [7, 8], "fulfil": [8, 11], "delta_": 8, "remap": 8, "onto": 8, "identifi": 8, "same": [7, 8, 10], "support": 8, "first": 8, "command": 8, "h1": [7, 8], "h1primalcel": [8, 10], "ne": 8, "take": [7, 8], "look": 8, "deform": [7, 8], "euler_angl": [3, 7, 8], "ones": [8, 10], "previou": 8, "e_k": 8, "e_0": 8, "e_1": 8, "canon": 8, "correct": 8, "covari": [7, 8], "transform": [7, 8, 10], "f_c": 8, "jacobian": [7, 8], "varphi": 8, "_c": 8, "definit": 8, "see": 8, "whole": [7, 8], "slightli": 8, "more": 8, "mesh_curl": 8, "435": 8, "hcurl": 8, "hcurldualcel": [7, 8, 10], "14": 8, "5714285714285716": 8, "gfu_curl": 8, "decompos": 8, "trace": [3, 8], "extens": 8, "minor": 8, "requir": [7, 8], "cube": [8, 10], "trilinear": 8, "semidiscret": 8, "written": 8, "down": 8, "similar": [7, 8, 10], "pleasant": 8, "check": 8, "By": 8, "face": [3, 8], "3p": 8, "12": 8, "6p": 8, "4p": 8, "verifi": [7, 8], "hold": 8, "unit_cub": [7, 8, 10], "hcurlprimalcel": [7, 8], "672": 8, "2736": 8, "7104": 8, "14640": 8, "count": 8, "nedg": 8, "nface": 8, "nelement": 8, "52": 8, "752": 8, "2964": 8, "7552": 8, "15380": 8, "would": 10, "chosen": 10, "classic": 10, "techniqu": 10, "x_p": 10, "yet": 10, "specifi": 10, "choos": 10, "exactli": 10, "quadratur": 10, "vanish": 10, "everi": 10, "except": 10, "expect": 10, "spars": 10, "gauss": 10, "radau": 10, "fix": 10, "treat": 10, "endpoint": 10, "gaussian": 10, "weight": 10, "easili": 10, "compar": 10, "exact": 10, "smooth": 10, "even": 10, "standard": 10, "fit": 10, "formula": 10, "trig_point": 10, "arrai": 10, "px": 10, "ob": 10, "spy": 10, "todens": 10, "seen": 10, "irs_f": 10, "m_diag": 10, "deletezeroel": [3, 10], "coupl": 10, "entri": 10, "fes_curl": 10, "irs_fes_curl": 10, "m_lump": 10, "store": [7, 10], "togeth": 10, "structur": 10, "less": 10, "exploit": 10, "access": 10, "m_exact": 10, "exacttim": 10, "m_superspars": 10, "stime": 10, "superspars": 10, "m_exact_inv": 10, "sparsecholeski": [3, 10], "m_supersparse_inv": 10, "tmp2": 10, "trialfunct": 7, "dx_vol": 7, "dx_edg": 7, "cross": [3, 7], "geometr": 7, "contribut": 7, "cancel": 7, "out": 7, "permut": 7, "equival": 7, "class": 7, "reduct": 7, "memori": 7, "cost": 7, "realiz": 7, "flag": 7, "curl_gf": 7, "curlt": 7, "gft": 7, "geometry_fre": 7, "100": [3, 7], "miss": 7, "spacial": 7, "nameerror": [], "befor": 7, "launch": 7, "hdiv": 7, "hdivdualcel": [], "h1primal": 7, "hdivdual": 7, "hiv": [], "side": 7, "mass_h1_inv": 7, "dx_h1": 7, "rh": 7, "gfp": 7, "1552574634552002": [], "014244794845581055": [], "16961169242858887": [], "029892921447753906": [], "15879392623901367": [], "013952970504760742": [], "15647482872009277": [], "024718046188354492": [], "attributeerror": [], "attribut": [], "15273332595825195": [], "014251708984375": [], "22423505783081055": [], "03263998031616211": [], "hprimaldualcel": [], "1410810947418213": [], "014180183410644531": [], "1572096347808838": [], "02831745147705078": [], "13854": [], "15579557418823242": [], "014397859573364258": [], "18164634704589844": [], "031702280044555664": [], "30510": [], "runtimeerror": [], "file": [], "401": [], "obj": [], "show": 3, "arg": [], "kwarg": [], "399": [], "402": [], "kwargs_with_default": [], "width": [], "height": [], "403": [], "404": [], "filenam": [], "405": [], "generatehtml": [], "widget": [], "76": [], "basewebguiscen": 3, "self": [], "74": [], "webguiwidget": [], "layout": [], "75": [], "encod": [], "getdata": [], "77": [], "displai": [], "260": [], "webglscen": [], "set_minmax": [], "257": [], "encodedata": [], "dtype": [], "float32": [], "259": [], "none": [], "rais": [], "cannot": [], "typ": [], "262": [], "263": [], "la": [], "basevector": [], "1414651870727539": [], "013903617858886719": [], "18274283409118652": [], "025954246520996094": [], "14151358604431152": [], "013034343719482422": [], "15631484985351562": [], "02606368064880371": [], "154870": 7, "344250": 7, "14334511756896973": [], "014346122741699219": [], "1761927604675293": [], "026534080505371094": [], "40": 7, "16558480262756348": [], "014064788818359375": [], "25002050399780273": [], "02897024154663086": [], "15956687927246094": [], "014201641082763672": [], "2367711067199707": [], "027669191360473633": [], "90": [], "13974499702453613": [], "01331949234008789": [], "1559600830078125": [], "024106740951538086": [], "150": 7, "14138245582580566": [], "013830423355102539": [], "20544767379760742": [], "028992652893066406": [], "1550123691558838": [], "014713525772094727": [], "22198820114135742": [], "026645660400390625": [], "includ": 7, "distribut": 7, "1411738395690918": [], "014023542404174805": [], "1531846523284912": [], "02953052520751953": [], "_h": 7, "int_t": 7, "mass_hdiv_inv": 7, "15952825546264648": [], "016004562377929688": [], "21116304397583008": [], "036969661712646484": [], "14396190643310547": [], "014764785766601562": [], "15488696098327637": [], "024342060089111328": [], "14777612686157227": [], "014744281768798828": [], "24978208541870117": [], "034796953201293945": [], "15024352073669434": [], "014367341995239258": [], "214324951171875": [], "03016209602355957": [], "m_u": [], "dot": [], "m_v": 11, "v_0": 11, "u_0": [], "u_": [], "u_j": [], "14087343215942383": [], "013223409652709961": [], "1841111183166504": [], "03317904472351074": [], "1401073932647705": [], "013926029205322266": [], "15143561363220215": [], "024031400680541992": [], "m_p": 11, "p_0": 11, "p_": 11, "p_j": 11, "13863372802734375": [], "014627456665039062": [], "2688016891479492": [], "02542257308959961": [], "known": 11, "sigma": 11, "spectral": 11, "radiu": [3, 11], "1457195281982422": [], "013231277465820312": [], "15278911590576172": [], "04138064384460449": [], "parital_t": [], "1678910255432129": [], "015158414840698242": [], "156968355178833": [], "029218196868896484": [], "14244389533996582": [], "013953447341918945": [], "2671644687652588": [], "03650474548339844": [], "14476919174194336": [], "013946533203125": [], "15114641189575195": [], "027591228485107422": [], "14116525650024414": [], "014034509658813477": [], "151580810546875": [], "02973484992980957": [], "locnr": [], "105662": [], "394338": [], "894338": [], "605662": [], "803561": [], "0982194": [], "0536948": [], "555556": [], "0778847": [], "0416667": [], "36656": [], "418661": [], "290669": [], "0296385": [], "725312": [], "0915628": [], "00999851": [], "51153": [], "0733767": [], "00616387": [], "341717": [], "391248": [], "0610607": [], "273846": [], "00365801": [], "316355": [], "227882": [], "0022025": [], "65232": [], "1840486526489258": [], "10226964950561523": [], "141755104064941": [], "0056302547454833984": [], "18722081184387207": [], "0008227825164794922": [], "14056038856506348": [], "014088630676269531": [], "15216398239135742": [], "030299901962280273": [], "1457281112670898": [], "0982968807220459": [], "94332766532898": [], "0051441192626953125": [], "1732780933380127": [], "000850677490234375": [], "341424e": [], "0x72bfd6786630": [], "0x72bfd4b248b0": [], "0x72bfd68f50b0": [], "44": [], "62706017494202": [], "707495e": [], "14597392082214355": [], "01627063751220703": [], "15703964233398438": [], "025795698165893555": [], "2047088146209717": [], "11965441703796387": [], "138391017913818": [], "005194664001464844": [], "17090559005737305": [], "0008075237274169922": [], "16479754447937012": [], "014780998229980469": [], "15702223777770996": [], "02642369270324707": [], "1339991092681885": [], "10081267356872559": [], "03842544555664": [], "0051000118255615234": [], "18155431747436523": [], "0008234977722167969": [], "1769404411315918": [], "020380258560180664": [], "17107176780700684": [], "024678468704223633": [], "122352": [], "913511037826538": [], "17615032196044922": [], "74958562850952": [], "008779048919677734": [], "4080026149749756": [], "003787517547607422": [], "14220356941223145": [], "013993024826049805": [], "2083144187927246": [], "03446626663208008": [], "0105292797088623": [], "21221518516540527": [], "43076729774475": [], "007737159729003906": [], "389972448348999": [], "0016677379608154297": [], "384611e": [], "0x7af0b864c4f0": [], "0x7af10877e570": [], "0x7af10877ecb0": [], "11096501350403": [], "603363e": [], "14337658882141113": [], "016633987426757812": [], "15613245964050293": [], "02317190170288086": [], "3169796466827393": [], "24944090843200684": [], "955610036849976": [], "007941484451293945": [], "39960145950317383": [], "0014188289642333984": [], "382613e": [], "0x7e7b5c13bcf0": [], "0x7e7b5dd9fbf0": [], "0x7e7b5dd9f0b0": [], "98267078399658": [], "421309e": [], "1393754482269287": [], "015882492065429688": [], "1510932445526123": [], "02279043197631836": [], "2654900550842285": [], "2561802864074707": [], "194032430648804": [], "008639097213745117": [], "3899693489074707": [], "0015819072723388672": [], "dualmesh": [0, 5], "gui": 3, "lx": 3, "ly": 3, "lz": 3, "obstacl": 3, "rad": 3, "experiment": 3, "t0": 3, "occ": 3, "wg": 3, "box": 3, "pnt": 3, "z": 3, "inflow": 3, "sphere": 3, "occgeometri": 3, "clip": 3, "120": 3, "15": 3, "fes_": 3, "hcurldualcells3d": 3, "fes_h": 3, "hcurlprimalcells3d": 3, "gf": 3, "dxe": 3, "dse": 3, "bf_mix": 3, "1280304": 3, "prepar": 3, "massh_inv": 3, "bfm_e": 3, "masse_inv": 3, "masse_surf": 3, "masse_surf_inv": 3, "getdof": 3, "inconsist": 3, "silenc": 3, "check_unus": 3, "plx": 3, "lfr": 3, "rhsefunc": 3, "lambda": 3, "01": 3, "scenebd": [], "draw_surf": 3, "sceneh": 3, "altshap": 3, "17777777777777762": [], "831430e": [], "35555555555555435": [], "598219e": [], "5333333333333311": [], "887032e": [], "7111111111111078": [], "012115e": [], "8888888888888845": [], "500": [], "951472e": [], "0666666666666653": [], "600": [], "115749e": [], "2444444444444531": [], "700": [], "067154e": [], "422222222222241": [], "800": [], "005547e": [], "6000000000000287": [], "900": [], "294784e": [], "7777777777778165": [], "093375e": [], "9555555555556043": [], "1100": [], "110797e": [], "133333333333392": [], "1200": [], "211048e": [], "3111111111111797": [], "1300": [], "056489e": [], "4888888888889675": [], "1400": [], "205722e": [], "6666666666667553": [], "1500": [], "146143e": [], "844444444444543": [], "1600": [], "149856e": [], "022222222222331": [], "1700": [], "331378e": [], "2000000000001187": [], "1800": [], "118499e": [], "3777777777779066": [], "1900": [], "036926e": [], "5555555555556944": [], "357804e": [], "733333333333482": [], "2100": [], "093240e": [], "91111111111127": [], "2200": [], "196624e": [], "19479799270629883": [], "02530384063720703": [], "2621617317199707": [], "06622910499572754": [], "gfe_anim": 3, "255827e": [], "ss": [], "13": [], "18": [], "19": [], "finsh": 3, "14386940002441406": [], "01591348648071289": [], "15303897857666016": [], "02896285057067871": [], "2036757469177246": [], "027491092681884766": [], "23320865631103516": [], "07455897331237793": [], "081141e": 3, "18926262855529785": 7, "02477884292602539": 7, "25863194465637207": 7, "08016109466552734": 7}, "objects": {}, "objtypes": {}, "objnames": {}, "titleterms": {"exampl": 0, "A": 1, "plane": 1, "wave": 1, "squar": 1, "ring": 2, "reson": [2, 9], "gener": [2, 8], "geometri": 2, "pml": 2, "sourc": 2, "coeffici": 2, "space": [2, 8, 9], "bilinear": 2, "form": 2, "time": [2, 6, 11], "loop": 2, "instal": 4, "get": 4, "ngsolv": [4, 5], "dual": [4, 5, 6, 8, 9], "cell": [4, 5, 6, 8, 9], "method": [3, 4, 5, 6, 9], "add": 4, "us": 4, "pip": 4, "work": 4, "code": 4, "test": 4, "dualcellspac": 4, "The": [5, 6, 7, 8], "tabl": 5, "content": 5, "refer": 5, "domain": 6, "maxwel": 6, "system": 6, "problem": [6, 9], "set": [6, 9], "basic": 9, "idea": 9, "state": 9, "art": 9, "from": 9, "fdtd": 9, "construct": [8, 9], "discret": [7, 8], "mass": 10, "lump": 10, "step": 11, "introduct": 9, "reason": 9, "explicit": 9, "scatter": 9, "type": 9, "galerkin": 9, "high": 9, "order": 9, "barycentr": 8, "mesh": 8, "refin": 8, "2d": 8, "map": 8, "basi": 8, "function": 8, "x_p": 8, "mathrm": 8, "grad": 8, "mathcal": 8, "t": 8, "tild": 8, "curl": 8, "summeri": 8, "three": 8, "dimens": 8, "assembl": [], "differenti": 7, "oper": 7, "etud": [], "comput": 7, "gradient": 7, "gaussian": 7, "peak": 7, "\u00e9tude": 7, "dualmesh": 3, "an": 3, "electromagnet": 3, "waveguid": 3, "3d": 3}, "envversion": {"sphinx.domains.c": 3, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 9, "sphinx.domains.index": 1, "sphinx.domains.javascript": 3, "sphinx.domains.math": 2, "sphinx.domains.python": 4, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx.ext.intersphinx": 1, "sphinxcontrib.bibtex": 9, "sphinx": 60}, "alltitles": {"Examples": [[0, "examples"]], "A plane wave on a square": [[1, "a-plane-wave-on-a-square"]], "Ring resonator": [[2, "ring-resonator"]], "Generating the geometry": [[2, "generating-the-geometry"]], "PML geometry": [[2, "pml-geometry"]], "Source": [[2, "source"]], "Coefficients": [[2, "coefficients"]], "Spaces": [[2, "spaces"]], "Bilinear forms": [[2, "bilinear-forms"]], "Time loop": [[2, "time-loop"]], "Dualmesh method for an electromagnetic waveguide in 3d": [[3, "dualmesh-method-for-an-electromagnetic-waveguide-in-3d"]], "Installation": [[4, "installation"]], "Get NGSolve": [[4, "get-ngsolve"]], "Install the Dual Cell Method add-on": [[4, "install-the-dual-cell-method-add-on"]], "Using pip": [[4, "using-pip"]], "Working with the code": [[4, "working-with-the-code"]], "Test the dualcellspaces installation": [[4, "test-the-dualcellspaces-installation"]], "The Dual Cell Method in NGSolve": [[5, "the-dual-cell-method-in-ngsolve"]], "Table of Contents": [[5, "table-of-contents"]], "References": [[5, "references"]], "The dual cell method for the time-domain Maxwell system": [[6, "the-dual-cell-method-for-the-time-domain-maxwell-system"]], "Problem setting": [[6, "problem-setting"]], "The discrete differential operators": [[7, "the-discrete-differential-operators"]], "\u00c9tude: computing the gradient of a Gaussian peak": [[7, "etude-computing-the-gradient-of-a-gaussian-peak"]], "Discrete spaces on dual cells": [[8, "discrete-spaces-on-dual-cells"]], "Barycentric mesh refinement in 2d": [[8, "barycentric-mesh-refinement-in-2d"]], "Mapped basis functions": [[8, "mapped-basis-functions"]], "The space X_P^{\\mathrm{grad}}(\\mathcal T)": [[8, "the-space-x-p-mathrm-grad-mathcal-t"]], "The space \\tilde X_P^{\\mathrm{curl}}(\\tilde{\\mathcal T})": [[8, "the-space-tilde-x-p-mathrm-curl-tilde-mathcal-t"]], "Summery of the construction": [[8, "summery-of-the-construction"]], "Generalization to three dimensions": [[8, "generalization-to-three-dimensions"]], "Introduction": [[9, "introduction"]], "Reasons for explicit methods": [[9, "reasons-for-explicit-methods"]], "scattering type problems": [[9, "scattering-type-problems"]], "resonance type problems": [[9, "resonance-type-problems"]], "State of the art": [[9, "state-of-the-art"]], "From FDTD to the dual cell method": [[9, "from-fdtd-to-the-dual-cell-method"]], "Basic idea of the dual cell construction": [[9, "basic-idea-of-the-dual-cell-construction"]], "Galerkin setting and high order spaces": [[9, "galerkin-setting-and-high-order-spaces"]], "Mass lumping": [[10, "mass-lumping"]], "Time stepping": [[11, "time-stepping"]]}, "indexentries": {}})