Search.setIndex({"docnames": ["examples", "examples/qualitative_dispersion", "examples/ring_resonator", "installation", "intro", "maxwell", "maxwell/introduction"], "filenames": ["examples.md", "examples/qualitative_dispersion.ipynb", "examples/ring_resonator.ipynb", "installation.md", "intro.md", "maxwell.md", "maxwell/introduction.md"], "titles": ["<span class=\"section-number\">3. </span>Examples", "<span class=\"section-number\">3.1. </span>A plane wave on a square", "<span class=\"section-number\">3.2. </span>Ring resonator", "<span class=\"section-number\">1. </span>Installation", "The Dual Cell Method in NGSolve", "<span class=\"section-number\">2. </span>The dual cell method for the time-domain Maxwell system", "<span class=\"section-number\">2.2. </span>Introduction"], "terms": {"A": [0, 4, 6], "plane": [0, 4], "wave": [0, 2, 4, 6], "squar": [0, 4], "ring": [0, 4], "reson": [0, 4], "we": [1, 2, 3, 4, 5, 6], "solv": [1, 2, 6], "two": [1, 2, 4, 6], "dimension": [1, 6], "equat": [1, 2, 4, 5], "find": [1, 5, 6], "h": [1, 2, 5, 6], "0": [1, 2, 5, 6], "t": [1, 2, 5, 6], "1": [1, 2, 6], "omega": [1, 2, 5, 6], "vector": [1, 2], "field": [1, 2, 4, 5, 6], "e": [1, 2, 3, 4, 5, 6], "mathrm": [1, 2, 5, 6], "div": [1, 2], "ar": [1, 2, 4, 5, 6], "begin": [1, 2, 5], "align": [1, 2, 5], "partial_t": [1, 2, 5, 6], "x": [1, 5], "nabla": [1, 2], "f": [1, 2], "exp": [1, 2], "400": [1, 2], "y": [1, 2], "2": [1, 2, 4, 6], "partial": [1, 6], "end": [1, 2, 5], "from": [1, 2, 3], "ngsolv": [1, 2], "import": [1, 2], "dualcellspac": [1, 2], "dc": [1, 2], "time": [1, 4, 6], "webgui": [1, 2], "draw": [1, 2], "after": 1, "necessari": [1, 6], "defin": [1, 2, 6], "some": [1, 5], "paramet": [1, 2, 6], "mesh": [1, 2, 4, 6], "maxh": [1, 2], "03": [1, 2], "tend": [1, 2], "order": [1, 2, 3, 4], "h0": 1, "cf": [1, 2], "20": [1, 2, 4], "e0": 1, "unit_squar": 1, "generatemesh": [1, 2], "space": 1, "fesh": 1, "h1dualcel": [1, 2], "fese": 1, "hdivprimalcel": [1, 2], "To": [1, 2, 5], "mass": [1, 2, 4, 6], "bilinear": 1, "form": [1, 6], "need": [1, 2, 5], "special": 1, "integr": [1, 2, 6], "rule": [1, 2], "de": 1, "tnt": [1, 2], "dh": 1, "dxh": 1, "dx": [1, 2], "intrul": [1, 2], "getintegrationrul": [1, 2], "dsw": [1, 2], "element_boundari": [1, 2], "true": [1, 2], "6": [1, 2], "dxw": [1, 2], "id": [1, 2], "massh": 1, "massinv": 1, "invers": [1, 2, 6], "massinvh": 1, "normal": [1, 2], "specialcf": [1, 2], "grad": 1, "bilinearform": [1, 2], "geom_fre": [1, 2], "assembl": [1, 2], "mat": [1, 2], "lffh": 1, "linearform": [1, 2], "The": [1, 2, 3, 6], "maxim": 1, "admiss": [1, 6], "step": [1, 6], "mai": [1, 2, 3, 5, 6], "estim": 1, "us": [1, 2, 6], "simpl": [1, 3], "power": 1, "iter": [1, 4], "def": [1, 2], "estimate_tau": 1, "maxstep": 1, "1000": 1, "tol": 1, "1e": 1, "4": [1, 2, 4], "vec": [1, 2], "createcolvector": 1, "setrandom": 1, "tmp": [1, 2], "createvector": [1, 2], "lam": 1, "i": [1, 2, 3, 4, 5, 6], "rang": 1, "print": [1, 2], "r": [1, 2, 5, 6], "data": [1, 2], "lamnew": 1, "innerproduct": 1, "tau": [1, 2, 6], "sqrt": [1, 2], "re": 1, "norm": 1, "diff": 1, "return": [1, 2], "did": 1, "converg": 1, "last": 1, "timestep": [1, 2, 6], "format": [1, 2], "9": [1, 2], "438193e": 1, "It": [1, 2, 4, 6], "remain": [1, 2, 6], "set": [1, 2], "initi": 1, "condit": [1, 5], "draweveri": [1, 2], "30": 2, "gfe": 1, "gridfunct": [1, 2], "gfh": 1, "gfh_histori": 1, "multidim": [1, 2], "scene": [1, 2], "intpoint": [1, 2], "getwebguipoint": [1, 2], "autoscal": [1, 2], "fals": [1, 2], "min": [1, 2], "max": [1, 2], "start": [1, 2, 6], "loop": 1, "now": [1, 2], "nowstart": 1, "energi": [1, 6], "tmph": 1, "tmpe": 1, "subtim": 1, "taskmanag": [1, 2], "while": [1, 2, 6], "timepass": 1, "before_energy_tim": 1, "addmultidimcompon": [1, 2], "redraw": [1, 2], "append": [1, 2], "current": 1, "dof": [1, 2], "": [1, 2, 6], "ndof": [1, 2], "comptim": 1, "n": [1, 2, 5], "per": 1, "second": 1, "1545": [], "1491700646156": [], "162256": 1, "653121e": [], "07": 2, "anim": [1, 2], "observ": 1, "preserv": 1, "modifi": [1, 3], "discret": [1, 2, 6], "matplotlib": [1, 2, 3], "pyplot": [1, 2], "pl": [1, 2], "plot": [1, 2], "ylim": 1, "simul": [2, 4], "propag": 2, "an": [2, 3, 4], "electromagnet": [2, 4], "tm": 2, "through": 2, "devic": 2, "dual": 2, "cell": [1, 2], "method": [1, 2], "add": 2, "desir": 2, "sketch": 2, "below": 2, "where": [2, 4, 5, 6], "horizont": 2, "waveguid": 2, "suppos": 2, "unbound": 2, "govern": 2, "scalar": [2, 6], "p": [2, 6], "u": 2, "l": [2, 4, 6], "int_": [2, 5], "cdot": [2, 6], "v": [2, 4, 6], "pq": 2, "q": [2, 6], "fq": 2, "all": [2, 5], "test": 2, "function": [1, 2, 4, 6], "domain": [2, 4, 6], "surround": 2, "absorb": 2, "layer": 2, "damp": 2, "system": [2, 4, 6], "alpha": 2, "omega_": 2, "tb": 2, "cup": 2, "lr": 2, "top": 2, "mathbf": [2, 5, 6], "hat": [2, 6], "left": [2, 6], "nn": 2, "right": [2, 6], "int": 2, "fv": 2, "omega_c": 2, "direct": 2, "case": 2, "given": [2, 6], "parallel": 2, "wire": 2, "shape": 2, "one": [2, 5, 6], "between": [2, 6], "splinegeometri": 2, "netgen": 2, "geom2d": 2, "gm": 2, "numpi": [2, 3], "np": 2, "geo": 2, "xneg": 2, "43": [2, 4], "xpo": 2, "yneg": 2, "48": 2, "ypo": 2, "wslab": 2, "04": 2, "cringx": 2, "cringi": 2, "rring": 2, "gap": 2, "005": 2, "pntx": 2, "pnty": 2, "pt": 2, "yi": 2, "xi": 2, "addpoint": 2, "inner": [2, 6], "rect": 2, "line": [2, 6], "leftdomain": 2, "rightdomain": 2, "3": [2, 5, 6], "5": [2, 4], "bc": 2, "normal_wg_rightbottom": 2, "normal_wg_leftbottom": 2, "7": 2, "normal_wg_righttop": 2, "8": 2, "normal_wg_lefttop": 2, "11": 2, "10": [2, 4], "addcircl": 2, "c": [2, 4, 5, 6], "setmateri": 2, "air": 2, "eps_nin": 2, "thi": [2, 4, 5, 6], "result": 2, "follow": [2, 6], "triangul": 2, "mesh_inn": 2, "curv": 2, "howev": [2, 6], "also": [2, 3, 5, 6], "want": [2, 3], "perfectli": 2, "match": 2, "pmlwidth": 2, "05": 2, "createpml": 2, "ha": [2, 3, 6], "ad": 2, "our": [2, 6], "interior": 2, "bd": 2, "getboundari": 2, "default": 2, "getmateri": 2, "pml_default": 2, "pml_corner": 2, "pml_default_duplicate_1": 2, "pml_normal_wg_rightbottom": 2, "pml_default_duplicate_2": 2, "pml_normal_wg_righttop": 2, "pml_default_duplicate_3": 2, "pml_default_duplicate_4": 2, "pml_default_duplicate_5": 2, "pml_normal_wg_lefttop": 2, "pml_default_duplicate_6": 2, "pml_normal_wg_leftbottom": 2, "pml_default_duplicate_7": 2, "boundari": [1, 2, 4, 5, 6], "wavelength": 2, "542": 2, "fcen": 2, "tpeak": 2, "sourcei": 2, "coefficientfunct": 2, "sin": 2, "pi": 2, "fes_facet": 2, "facetfespac": 2, "gfsourc": 2, "definedon": 2, "t_envelop": 2, "ab": 2, "els": 2, "delta": 2, "001": 2, "arang": 2, "xlabel": 2, "figur": 2, "materi": 2, "background": 2, "medium": [2, 5], "eps_r": 2, "startswith": 2, "pml_normal_wg": 2, "ep": 2, "scale": [2, 6], "non": 2, "corner": 2, "nvec": 2, "cfn": 2, "next": 2, "up": 2, "fes_u": 2, "fes_p": 2, "dirichlet": 2, "outer": 2, "vectori": [2, 6], "fe": 2, "total": 2, "030700e": [], "887000e": [], "383540e": [], "which": [2, 6], "respect": [2, 4, 6], "These": 2, "can": [2, 3], "obtain": [2, 6], "via": [2, 3, 4], "over": [1, 2], "volum": 2, "element": [1, 2, 3, 4, 6], "ir": 2, "et": 2, "segm": 2, "fem": 2, "integrationrul": 2, "object": 2, "0x70232cb071f0": [], "trig": 2, "0x70232cb06930": [], "tet": 2, "0x70232cb07630": [], "diverg": 2, "contain": 2, "jump": 2, "term": [1, 2, 6], "due": [2, 6], "fact": [2, 6], "v_h": 2, "discontinu": [2, 4, 6], "b": [2, 4, 5], "matric": [2, 4], "block": [1, 2, 4], "diagon": [1, 2, 4], "lump": [2, 4], "thei": 2, "implement": [2, 4, 6], "effici": 2, "fespac": 2, "invmassp": 2, "freedof": 2, "invmassu": 2, "dim": 2, "pml1d": 2, "pml_normal": 2, "dampp1": 2, "dampp2": 2, "dampu1": 2, "outerproduct": 2, "dampu2": 2, "big": 2, "oper": [2, 6], "embed": 2, "small": 2, "emb_p": 2, "emb_phat": 2, "emb_u": 2, "emb_uhat": 2, "b_big": 2, "dampu_big": 2, "dampp_big": 2, "invmassp_big": 2, "invmassu_big": 2, "lastli": 2, "q_big": [], "testfunct": [], "lsrc": 2, "gfu": 2, "gfp_histori": 2, "compon": [2, 6], "4e": [], "5e": [], "25": [], "startnow": 2, "200": [], "mdofss": [], "087120000001544": [], "113300e": [], "ssss": [], "keyboardinterrupt": [], "traceback": [], "most": [3, 6], "recent": [3, 6], "call": [], "In": 5, "37": [], "23": 4, "21": [], "22": [], "24": [], "finish": [], "packag": 3, "high": [3, 4], "finit": [3, 4, 6], "librari": [3, 4], "thu": [3, 4, 6], "main": [3, 6], "premis": 3, "suffici": 3, "version": 3, "beforehand": 3, "wai": 3, "python": 3, "m": [3, 4], "scipi": 3, "jupyt": 3, "ipyparallel": 3, "scikit": 3, "build": 3, "upgrad": 3, "webgui_jupyter_widget": 3, "For": [3, 4, 6], "troubleshoot": 3, "refer": 3, "variou": 3, "tutori": 3, "ngs24": 3, "document": 3, "If": 3, "you": 3, "do": 3, "updat": 3, "your": 3, "might": 3, "consid": 3, "virtual": 3, "environ": 3, "As": 3, "avail": 3, "core": 3, "again": [3, 6], "have": [3, 6], "binari": 3, "git": 3, "http": 3, "github": 3, "com": 3, "dcm": [3, 4], "built": 3, "sourc": [1, 3, 5], "dev": 3, "wa": 3, "g": [3, 4, 6], "pre": 3, "pybind11_stubgen": 3, "isol": 3, "modul": 3, "clone": 3, "them": 3, "done": [3, 6], "either": 3, "cmake": 3, "cd": 3, "mkdir": 3, "make": [3, 6], "j4": 3, "whether": 3, "demo": 3, "dc_intrul": 3, "test_spac": 3, "wess": 4, "j": [4, 6], "sch\u00f6berl": 4, "tu": 4, "wien": 4, "institut": 4, "analysi": 4, "scientif": 4, "comput": [2, 4], "base": 4, "joint": 4, "work": 4, "kapidani": 4, "codecasa": 4, "book": 4, "design": 4, "provid": 4, "introduct": 4, "exampl": 4, "galerkin": 4, "acoustiv": 4, "mix": 4, "formul": [4, 5, 6], "disconitinu": 4, "variant": 4, "approxim": [4, 6], "conform": 4, "each": [4, 6], "other": [4, 6], "ansatz": 4, "featur": [4, 6], "differ": [4, 6], "full": 4, "mathemat": 4, "kcschoberl21": 4, "wkcs23": 4, "instal": 4, "maxwel": 4, "bernard": 4, "lorenzo": 4, "joachim": 4, "sch": 4, "\u00f6": 4, "berl": 4, "arbitrari": 4, "depend": [4, 6], "d": [4, 5], "axwel": 4, "phy": 4, "433": 4, "paper": 4, "No": 4, "110184": 4, "2021": 4, "doi": 4, "1016": 4, "jcp": 4, "marku": 4, "polynomi": 4, "degre": [4, 6], "acoust": 4, "2023": 4, "arxiv": 4, "2312": 4, "14716": 4, "section": [5, 6], "describ": 5, "detail": 5, "how": 5, "appli": [5, 6], "neglect": 5, "possibl": [5, 6], "here": 5, "curl": 5, "text": [5, 6], "suitabl": [1, 5], "subset": 5, "mathbb": [5, 6], "close": 5, "constitut": 5, "relat": 5, "color": 5, "emph1": 5, "varepsilon": [5, 6], "emph2": 5, "mu": [5, 6], "permitt": 5, "permeabl": 5, "question": 5, "weak": [5, 6], "eh": 5, "assum": 5, "homogen": 5, "abov": 5, "rewritten": 5, "377423e": [], "08": 1, "270705e": [], "187766e": [], "version_info": [], "modulenotfounderror": [], "input": [], "name": [], "9933350511855035": [], "1540": [], "06247129244141534": [], "381633e": [], "088": [], "160984e": [], "392659e": [], "35": [], "409019e": [], "322534e": [], "371604e": [], "0x7abf0dc794b0": [], "0x7abf0dc793f0": [], "0x7abf0dc78df0": [], "036663e": [], "3e": 2, "300": [], "dofss": 2, "when": 6, "type": [], "problem": [], "approach": 6, "typic": 6, "choic": 6, "implicit": 6, "explicit": [], "former": 6, "matrix": 6, "compos": 6, "differenti": [1, 6], "although": 6, "scheme": 6, "uncondition": 6, "stabl": 6, "independ": 6, "factor": 6, "larg": 6, "doe": 6, "well": 6, "veri": 6, "number": 6, "freedom": 6, "hand": 6, "mere": 6, "There": 6, "exist": 6, "sever": 6, "conveni": 6, "downsid": 6, "stabil": 6, "qualiti": 6, "finer": 6, "largest": 6, "guarante": 6, "popular": 6, "afdsj": [], "adf": [], "dg": 6, "appelogr20": [4, 6], "daniel": 4, "appel": 4, "fortino": 4, "garcia": 4, "olof": 4, "runborg": 4, "waveholtz": 4, "solut": 4, "helmholtz": [4, 6], "siam": 4, "journal": 4, "42": 4, "a1950": 4, "a1983": 4, "2020": [4, 6], "1137": 4, "19m1299062": 4, "apart": 6, "obviou": 6, "better": 6, "applic": 6, "solver": [4, 6], "interest": 6, "frequenc": 6, "scatter": [], "precondition": [4, 6], "stolk": 4, "eigenvalu": [4, 6], "lite": 6, "pinvit": 6, "lobpcg": 6, "earli": 6, "2000": 6, "knyazev": 4, "filter": [4, 6], "2e5": 6, "hit": 6, "googl": 6, "scholar": 6, "2e4": 6, "sinc": 6, "2022": 6, "low": 6, "onli": 6, "hexahedr": 6, "grid": 6, "5e4": 6, "6e3": 6, "numer": 6, "flux": 6, "penalti": 6, "both": 6, "been": 6, "around": 6, "late": 6, "1960": 6, "1970": 6, "still": 6, "wide": 6, "engin": 6, "sto21": [4, 6], "christiaan": 4, "a3469": 4, "a3502": 4, "20m1359997": 4, "kny01": [4, 6], "andrew": 4, "toward": 4, "optim": [1, 4], "precondit": 4, "eigensolv": 4, "local": [4, 6], "conjug": 4, "gradient": 4, "517": 4, "541": 4, "2001": 4, "s1064827500366124": 4, "nw24": [4, 6], "lothar": 4, "nannen": 4, "krylov": 4, "2024": 4, "2402": 4, "08515": 4, "idea": [], "stoke": 6, "theorem": 6, "my": [], "caption": [], "visual": 6, "lead": 6, "quantiti": 6, "interlac": 6, "pointwis": 6, "mimmick": 6, "l_": 6, "k": 6, "point": [1, 6], "valu": 6, "x_i": 6, "y_j": 6, "z_k": 6, "tangenti": 6, "stagger": 6, "h_": 6, "tfrac": 6, "approx": 6, "e_": 6, "satisfi": 6, "leap": 6, "frog": 6, "size": 6, "_": 6, "frac": 6, "quad": 6, "carri": 6, "goal": 6, "gener": 6, "tetrahedr": 6, "157738e": [], "basi": [1, 6], "focuss": 6, "go": 6, "dimens": 6, "primal": [1, 6], "interpret": 6, "orang": 6, "circl": 6, "piecewis": 6, "constant": 6, "dark": 6, "grei": 6, "consist": 6, "four": 6, "proce": 6, "unknown": 6, "neighbour": 6, "continu": 6, "tild": 6, "c_j": 6, "bigcup_": 6, "c_": 6, "mathcal": 6, "_j": 6, "w_": 6, "setr": [], "foral": 6, "cont": 6, "v_": 6, "global": 6, "w": 6, "zero": 6, "outsid": 6, "origin": 6, "pose": 6, "sum": [], "rot": 6, "sum_": 6, "semi": 6, "ultra": 6, "variat": 6, "remark": 6, "aris": 6, "across": 6, "part": 6, "help": 6, "skew": 6, "symmetri": 6, "immedi": 6, "conserv": 6, "simplex": 6, "motiv": 6, "construct": [], "symbol": 1, "correspond": 1, "nodal": 1, "higher": 1, "sparsiti": 1, "pattern": 1, "328608e": [], "596200e": 2, "352700e": 2, "824640e": 2, "0x7efeaab2eef0": [], "0x7efeaab2f5f0": [], "0x7efeaab2c9f0": [], "879820000000851": [], "030279e": [], "sss": [], "averag": 2, "1546": 1, "424357e": [], "home": [], "mwess": [], "lib": [], "python3": [], "site": [], "project": [], "__init__": [], "py": [], "63": [], "userwarn": [], "unabl": [], "axes3d": [], "multipl": [], "being": [], "pip": [], "3d": [], "warn": [], "0x737388a1c9b0": [], "0x737388a1c270": [], "0x73738c7f48b0": [], "11327": 2, "613625e": [], "309617e": [], "0x7c895f45ba70": [], "0x7c895d73d2b0": [], "0x7c895f343330": [], "552361e": [], "50": 1, "342434e": [], "498972e": [], "265443e": 1, "0x7f00a4f15030": 2, "0x7f00a4fa5530": 2, "0x7f00a4fa5d30": 2, "735248e": 2}, "objects": {}, "objtypes": {}, "objnames": {}, "titleterms": {"exampl": 0, "A": 1, "plane": 1, "wave": 1, "squar": 1, "ring": 2, "reson": [2, 6], "gener": 2, "geometri": 2, "pml": 2, "sourc": 2, "coeffici": 2, "space": [2, 6], "bilinear": 2, "form": 2, "time": [2, 5], "loop": 2, "instal": 3, "get": 3, "ngsolv": [3, 4], "dual": [3, 4, 5, 6], "cell": [3, 4, 5, 6], "method": [3, 4, 5, 6], "add": 3, "us": 3, "pip": 3, "work": 3, "code": 3, "test": 3, "dualcellspac": 3, "The": [4, 5], "tabl": 4, "content": 4, "refer": 4, "domain": 5, "maxwel": 5, "system": 5, "problem": [5, 6], "set": [5, 6], "basic": 6, "idea": 6, "state": 6, "art": 6, "from": 6, "fdtd": 6, "construct": 6, "discret": [], "mass": [], "lump": [], "step": [], "introduct": 6, "reason": 6, "explicit": 6, "scatter": 6, "type": 6, "galerkin": 6, "high": 6, "order": 6}, "envversion": {"sphinx.domains.c": 3, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 9, "sphinx.domains.index": 1, "sphinx.domains.javascript": 3, "sphinx.domains.math": 2, "sphinx.domains.python": 4, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx.ext.intersphinx": 1, "sphinxcontrib.bibtex": 9, "sphinx": 60}, "alltitles": {"Examples": [[0, "examples"]], "A plane wave on a square": [[1, "a-plane-wave-on-a-square"]], "Ring resonator": [[2, "ring-resonator"]], "Generating the geometry": [[2, "generating-the-geometry"]], "PML geometry": [[2, "pml-geometry"]], "Source": [[2, "source"]], "Coefficients": [[2, "coefficients"]], "Spaces": [[2, "spaces"]], "Bilinear forms": [[2, "bilinear-forms"]], "Time loop": [[2, "time-loop"]], "Installation": [[3, "installation"]], "Get NGSolve": [[3, "get-ngsolve"]], "Install the Dual Cell Method add-on": [[3, "install-the-dual-cell-method-add-on"]], "Using pip": [[3, "using-pip"]], "Working with the code": [[3, "working-with-the-code"]], "Test the dualcellspaces installation": [[3, "test-the-dualcellspaces-installation"]], "The Dual Cell Method in NGSolve": [[4, "the-dual-cell-method-in-ngsolve"]], "Table of Contents": [[4, "table-of-contents"]], "References": [[4, "references"]], "The dual cell method for the time-domain Maxwell system": [[5, "the-dual-cell-method-for-the-time-domain-maxwell-system"]], "Problem setting": [[5, "problem-setting"]], "Introduction": [[6, "introduction"]], "Reasons for explicit methods": [[6, "reasons-for-explicit-methods"]], "scattering type problems": [[6, "scattering-type-problems"]], "resonance type problems": [[6, "resonance-type-problems"]], "State of the art": [[6, "state-of-the-art"]], "From FDTD to the dual cell method": [[6, "from-fdtd-to-the-dual-cell-method"]], "Basic idea of the dual cell construction": [[6, "basic-idea-of-the-dual-cell-construction"]], "Galerkin setting and high order spaces": [[6, "galerkin-setting-and-high-order-spaces"]]}, "indexentries": {}})