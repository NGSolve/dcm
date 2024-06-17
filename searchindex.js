Search.setIndex({"docnames": ["examples", "examples/qualitative_dispersion", "examples/ring_resonator", "examples/waveguide_3d", "installation", "intro", "maxwell", "maxwell/bibliography", "maxwell/diffops", "maxwell/hcurl", "maxwell/introduction", "maxwell/mass_lumping", "maxwell/time_stepping"], "filenames": ["examples.md", "examples/qualitative_dispersion.ipynb", "examples/ring_resonator.ipynb", "examples/waveguide_3d.ipynb", "installation.md", "intro.md", "maxwell.md", "maxwell/bibliography.md", "maxwell/diffops.md", "maxwell/hcurl.md", "maxwell/introduction.md", "maxwell/mass_lumping.md", "maxwell/time_stepping.md"], "titles": ["<span class=\"section-number\">3. </span>Examples", "<span class=\"section-number\">3.1. </span>A plane wave on a square", "<span class=\"section-number\">3.2. </span>Ring resonator", "<span class=\"section-number\">3.3. </span>Electromagnetic 3d waveguide", "<span class=\"section-number\">1. </span>Installation", "The Dual Cell Method in NGSolve", "<span class=\"section-number\">2. </span>The dual cell method for the time-domain Maxwell system", "<span class=\"section-number\">2.7. </span>Bibliography", "<span class=\"section-number\">2.5. </span>The discrete differential operators", "<span class=\"section-number\">2.3. </span>Discrete spaces on dual cells", "<span class=\"section-number\">2.2. </span>Introduction", "<span class=\"section-number\">2.4. </span>Mass lumping", "<span class=\"section-number\">2.6. </span>Time stepping"], "terms": {"A": [0, 5, 7, 9, 10], "plane": [0, 3, 5], "wave": [0, 2, 5, 7, 8, 10], "squar": [0, 5, 9, 11], "ring": [0, 5], "reson": [0, 5], "we": [1, 2, 4, 5, 6, 8, 9, 10, 11, 12], "solv": [1, 2, 8, 10], "two": [1, 2, 5, 9, 10], "dimension": [1, 9, 10], "equat": [1, 2, 5, 6, 7, 12], "find": [1, 6, 8, 9, 10], "h": [1, 2, 3, 6, 8, 9, 10], "0": [1, 2, 3, 6, 8, 9, 10, 11, 12], "t": [1, 2, 3, 6, 8, 10], "1": [1, 2, 3, 8, 9, 10, 11, 12], "omega": [1, 2, 3, 6, 9, 10], "vector": [1, 2, 8, 9], "field": [1, 2, 5, 6, 10], "e": [1, 2, 3, 4, 5, 6, 8, 9, 10], "mathrm": [1, 2, 6, 8, 10], "div": [1, 2, 8], "ar": [1, 2, 5, 6, 8, 9, 10, 11, 12], "begin": [1, 2, 6, 12], "align": [1, 2, 6, 12], "partial_t": [1, 2, 6, 9, 10, 12], "x": [1, 3, 6, 8, 9], "nabla": [1, 2], "f": [1, 2, 8, 9], "exp": [1, 2, 3, 8], "400": [1, 2], "y": [1, 2, 3, 8], "2": [1, 2, 3, 5, 7, 8, 9, 10, 11, 12], "partial": [1, 8, 9, 10], "end": [1, 2, 3, 6, 9, 12], "from": [1, 2, 3, 4, 8, 9, 11], "ngsolv": [1, 2, 3, 8, 9, 11], "import": [1, 2, 3, 8, 9, 11], "dualcellspac": [1, 2, 3, 8, 9, 11], "dc": [1, 2, 3, 8, 9, 11], "time": [1, 3, 5, 7, 8, 9, 10, 11], "webgui": [1, 2, 3, 8, 9], "draw": [1, 2, 3, 8, 9], "after": 1, "necessari": [1, 10], "defin": [1, 2, 8, 9, 10, 11], "some": [1, 6, 8], "paramet": [1, 2, 3, 10], "mesh": [1, 2, 3, 5, 8, 10, 11], "maxh": [1, 2, 3, 8, 9, 11], "03": [1, 2, 8], "tend": [1, 2, 3], "order": [1, 2, 3, 4, 5, 7, 8, 9, 11], "h0": 1, "cf": [1, 2, 3, 8], "20": [1, 2, 5, 7], "e0": 1, "unit_squar": [1, 8, 9, 11], "generatemesh": [1, 2, 3, 8, 9, 11], "space": [1, 8, 11], "fesh": [1, 8], "h1dualcel": [1, 2, 8], "fese": [1, 8], "hdivprimalcel": [1, 2, 8], "To": [1, 2, 6, 9, 11], "mass": [1, 2, 3, 5, 7, 8, 10, 12], "bilinear": [1, 3, 9], "form": [1, 3, 10, 12], "need": [1, 2, 6, 8, 11], "special": [1, 3], "integr": [1, 2, 3, 8, 9, 10, 11], "rule": [1, 2, 3, 8, 11], "de": [1, 3, 8], "tnt": [1, 2, 3, 8, 11], "dh": [1, 3], "dxh": [1, 3], "dx": [1, 2, 3, 8, 11], "intrul": [1, 2, 3, 8, 11], "getintegrationrul": [1, 2, 3, 8, 11], "dsw": [1, 2, 3, 8], "element_boundari": [1, 2, 3, 8], "true": [1, 2, 3, 8, 9], "6": [1, 2, 9], "dxw": [1, 2, 3, 8], "id": [1, 2], "massh": 1, "massinv": 1, "invers": [1, 2, 3, 8, 10, 11], "massinvh": 1, "normal": [1, 2, 3, 8], "specialcf": [1, 2, 3, 8], "grad": [1, 8], "bilinearform": [1, 2, 3, 8, 11], "geom_fre": [1, 2, 3, 8], "assembl": [1, 2, 3, 8, 11], "mat": [1, 2, 3, 8, 11], "lffh": 1, "linearform": [1, 2, 3, 8], "The": [1, 2, 4, 10, 11], "maxim": 1, "admiss": [1, 10], "step": [1, 3, 9, 10, 11], "mai": [1, 2, 4, 6, 9, 10, 11], "estim": 1, "us": [1, 2, 3, 8, 9, 10, 11, 12], "simpl": [1, 4, 9], "power": 1, "iter": [1, 7], "def": [1, 2], "estimate_tau": 1, "maxstep": 1, "1000": 1, "tol": 1, "1e": [1, 3, 11], "4": [1, 2, 3, 7, 8, 9, 12], "vec": [1, 2, 3, 8, 9], "createcolvector": 1, "setrandom": [1, 8, 11], "tmp": [1, 2, 11], "createvector": [1, 2, 3, 11], "lam": 1, "i": [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12], "rang": [1, 8, 9, 11], "print": [1, 2, 3, 8, 9, 11], "r": [1, 2, 3, 6, 8, 9, 10], "data": [1, 2, 3, 8, 11], "lamnew": 1, "innerproduct": [1, 3], "tau": [1, 2, 3, 10, 12], "sqrt": [1, 2], "re": 1, "norm": 1, "diff": 1, "return": [1, 2], "did": 1, "converg": 1, "last": 1, "timestep": [1, 2, 3, 10, 12], "format": [1, 2, 3, 8, 9, 11], "9": [1, 2], "438193e": 1, "It": [1, 2, 5, 10, 12], "remain": [1, 2, 10], "set": [1, 2, 3, 9], "initi": 1, "condit": [1, 6], "draweveri": [1, 2, 3], "30": 2, "gfe": [1, 3, 8], "gridfunct": [1, 2, 3, 8, 9], "gfh": [1, 3, 8], "gfh_histori": 1, "multidim": [1, 2, 3, 9], "scene": [1, 2, 3], "intpoint": [1, 2, 9], "getwebguipoint": [1, 2, 3, 8, 9], "autoscal": [1, 2, 3], "fals": [1, 2, 3], "min": [1, 2, 3, 9], "max": [1, 2, 3, 9], "start": [1, 2, 9, 10, 11], "loop": 1, "now": [1, 2, 3, 8, 11], "nowstart": 1, "energi": [1, 10], "tmph": 1, "tmpe": [1, 3], "subtim": 1, "taskmanag": [1, 2, 3, 8, 11], "while": [1, 2, 3, 10], "timepass": [1, 3], "before_energy_tim": 1, "addmultidimcompon": [1, 2, 3], "redraw": [1, 2, 3], "append": [1, 2], "current": 1, "dof": [1, 2, 3, 8, 9], "": [1, 2, 3, 8, 9, 10, 11], "ndof": [1, 2, 3, 8, 9, 11], "comptim": 1, "n": [1, 2, 6, 8, 9, 11], "per": [1, 9], "second": [1, 2], "1545": 1, "1491700646156": 1, "162256": 1, "653121e": [], "07": [2, 3], "anim": [1, 2, 3, 9], "observ": [1, 11], "preserv": 1, "modifi": [1, 4], "discret": [1, 2, 10, 12], "matplotlib": [1, 2, 3, 4, 11], "pyplot": [1, 2, 3, 11], "pl": [1, 2, 3, 11], "plot": [1, 2, 3, 11], "ylim": 1, "simul": [2, 5], "propag": 2, "an": [2, 3, 4, 5, 7, 9], "electromagnet": [0, 2, 5, 7], "tm": 2, "through": 2, "devic": 2, "dual": [2, 7, 11], "cell": [1, 2, 7, 8, 11], "method": [1, 2, 7, 11], "add": 2, "desir": 2, "sketch": 2, "below": 2, "where": [2, 5, 6, 9, 10, 11, 12], "horizont": 2, "waveguid": [0, 2, 5], "suppos": 2, "unbound": 2, "govern": 2, "scalar": [2, 10], "p": [2, 8, 9, 10, 11, 12], "u": [2, 8, 11], "l": [2, 5, 9, 10], "int_": [2, 6, 8], "cdot": [2, 8, 9, 10], "v": [2, 7, 8, 10, 11, 12], "pq": 2, "q": [2, 8, 9, 10], "fq": 2, "all": [2, 6, 8, 9], "test": 2, "function": [1, 2, 5, 8, 10, 11], "domain": [2, 5, 7, 8, 9, 10], "surround": 2, "absorb": 2, "layer": 2, "damp": 2, "system": [2, 5, 10, 12], "alpha": 2, "omega_": 2, "tb": 2, "cup": 2, "lr": 2, "top": [2, 8, 9, 12], "mathbf": [2, 6, 8, 9, 10, 12], "hat": [2, 9, 10], "left": [2, 9, 10, 12], "nn": 2, "right": [2, 8, 9, 10, 12], "int": [2, 9], "fv": 2, "omega_c": 2, "direct": [2, 9], "case": 2, "given": [2, 9, 10, 12], "parallel": 2, "wire": 2, "shape": [2, 11], "one": [2, 6, 8, 9, 10, 11], "between": [2, 10, 11], "splinegeometri": 2, "netgen": [2, 3], "geom2d": 2, "gm": 2, "numpi": [2, 3, 4, 11], "np": [2, 3, 11], "geo": [2, 3], "xneg": 2, "43": [2, 7], "xpo": 2, "yneg": 2, "48": [2, 9], "ypo": 2, "wslab": 2, "04": 2, "cringx": 2, "cringi": 2, "rring": 2, "gap": 2, "005": 2, "pntx": 2, "pnty": 2, "pt": 2, "yi": 2, "xi": [2, 9], "addpoint": 2, "inner": [2, 9, 10], "rect": 2, "line": [2, 10], "leftdomain": 2, "rightdomain": 2, "3": [2, 3, 6, 8, 9, 10, 11], "5": [2, 3, 7, 8, 9, 11], "bc": 2, "normal_wg_rightbottom": 2, "normal_wg_leftbottom": 2, "7": [2, 9], "normal_wg_righttop": 2, "8": [2, 3], "normal_wg_lefttop": 2, "11": 2, "10": [2, 3, 5, 7, 11], "addcircl": 2, "c": [2, 3, 6, 7, 9, 10], "setmateri": 2, "air": 2, "eps_nin": 2, "thi": [2, 3, 5, 6, 8, 9, 10, 11, 12], "result": [2, 8, 11], "follow": [2, 8, 9, 10], "triangul": 2, "mesh_inn": 2, "curv": [2, 3], "howev": [2, 10], "also": [2, 4, 6, 10], "want": [2, 4], "perfectli": 2, "match": 2, "pmlwidth": 2, "05": 2, "createpml": 2, "ha": [2, 4, 9, 10], "ad": 2, "our": [2, 8, 9, 10, 11], "interior": 2, "bd": 2, "getboundari": 2, "default": 2, "getmateri": 2, "pml_default": 2, "pml_corner": 2, "pml_default_duplicate_1": 2, "pml_normal_wg_rightbottom": 2, "pml_default_duplicate_2": 2, "pml_normal_wg_righttop": 2, "pml_default_duplicate_3": 2, "pml_default_duplicate_4": 2, "pml_default_duplicate_5": 2, "pml_normal_wg_lefttop": 2, "pml_default_duplicate_6": 2, "pml_normal_wg_leftbottom": 2, "pml_default_duplicate_7": 2, "boundari": [1, 2, 3, 5, 6, 8, 9, 10], "wavelength": 2, "542": 2, "fcen": 2, "tpeak": 2, "sourcei": 2, "coefficientfunct": 2, "sin": [2, 3], "pi": [2, 3], "fes_facet": 2, "facetfespac": 2, "gfsourc": 2, "definedon": [2, 3], "t_envelop": 2, "ab": 2, "els": 2, "delta": 2, "001": 2, "arang": [2, 3], "xlabel": 2, "figur": [2, 11], "materi": 2, "background": 2, "medium": [2, 6], "eps_r": 2, "startswith": 2, "pml_normal_wg": 2, "ep": 2, "scale": [2, 10], "non": 2, "corner": 2, "nvec": 2, "cfn": 2, "next": [2, 9], "up": [2, 8, 9], "fes_u": 2, "fes_p": 2, "dirichlet": [2, 3], "outer": 2, "vectori": [2, 9, 10, 11], "fe": [2, 3, 11], "total": [2, 3], "030700e": [], "887000e": [], "383540e": [], "which": [2, 8, 9, 10, 11], "respect": [2, 5, 9, 10, 11], "These": 2, "can": [2, 4, 9, 11], "obtain": [2, 8, 9, 10, 11], "via": [2, 4, 7, 9, 11], "over": [1, 2, 3, 8], "volum": 2, "element": [1, 2, 4, 5, 8, 9, 10, 11], "ir": [2, 8, 11], "et": [2, 11], "segm": 2, "fem": 2, "integrationrul": [2, 11], "object": [2, 3], "0x70232cb071f0": [], "trig": [2, 11], "0x70232cb06930": [], "tet": 2, "0x70232cb07630": [], "diverg": 2, "contain": [2, 9], "jump": 2, "term": [1, 2, 8, 10], "due": [2, 8, 10], "fact": [2, 10, 11], "v_h": 2, "discontinu": [2, 5, 10, 11], "b": [2, 5, 6, 8, 12], "matric": [2, 5, 7, 8, 11, 12], "block": [1, 2, 5, 7, 11], "diagon": [1, 2, 5, 7, 11], "lump": [2, 5, 7, 8], "thei": [2, 11], "implement": [2, 5, 8, 9, 10, 11], "effici": 2, "fespac": [2, 11], "invmassp": 2, "freedof": [2, 3], "invmassu": 2, "dim": 2, "pml1d": 2, "pml_normal": 2, "dampp1": 2, "dampp2": 2, "dampu1": 2, "outerproduct": 2, "dampu2": 2, "big": 2, "oper": [2, 3, 10, 12], "embed": 2, "small": 2, "emb_p": 2, "emb_phat": 2, "emb_u": 2, "emb_uhat": 2, "b_big": 2, "dampu_big": 2, "dampp_big": 2, "invmassp_big": 2, "invmassu_big": 2, "lastli": [2, 8], "q_big": [], "testfunct": 8, "lsrc": 2, "gfu": [2, 8, 9], "gfp_histori": 2, "compon": [2, 3, 10], "4e": [], "5e": [], "25": 9, "startnow": 2, "200": [], "mdofss": [], "087120000001544": [], "113300e": [], "ssss": [], "keyboardinterrupt": [], "traceback": [], "most": [4, 10], "recent": [4, 10], "call": 9, "In": [6, 9, 11], "37": [], "23": 7, "21": [], "22": [], "24": [], "finish": [], "packag": [4, 9, 11], "high": [4, 5, 11], "finit": [4, 5, 9, 10, 11], "librari": [3, 4, 5], "thu": [4, 5, 8, 9, 10, 11], "main": [4, 10, 11], "premis": 4, "suffici": 4, "version": 4, "beforehand": 4, "wai": [4, 8, 9], "python": 4, "m": [4, 5, 7, 11], "scipi": 4, "jupyt": 4, "ipyparallel": 4, "scikit": 4, "build": 4, "upgrad": 4, "webgui_jupyter_widget": 4, "For": [4, 5, 9, 10, 11, 12], "troubleshoot": 4, "refer": [4, 8, 9, 11], "variou": 4, "tutori": 4, "ngs24": 4, "document": 4, "If": 4, "you": 4, "do": 4, "updat": 4, "your": 4, "might": 4, "consid": 4, "virtual": 4, "environ": 4, "As": [4, 8, 11], "avail": 4, "core": 4, "again": [4, 9, 10], "have": [4, 9, 10, 11], "binari": 4, "git": 4, "http": [4, 5, 7], "github": 4, "com": 4, "dcm": [4, 5], "built": 4, "sourc": [1, 4, 6], "dev": 4, "wa": [4, 11], "g": [4, 5, 10], "pre": 4, "pybind11_stubgen": 4, "isol": 4, "modul": 4, "clone": 4, "them": [4, 9], "done": [4, 8, 9, 10], "either": 4, "cmake": 4, "cd": 4, "mkdir": 4, "make": [4, 10], "j4": 4, "whether": [4, 11], "demo": 4, "dc_intrul": 4, "test_spac": 4, "wess": [5, 7], "j": [5, 7, 9, 10, 12], "sch\u00f6berl": [5, 7], "tu": 5, "wien": 5, "institut": 5, "analysi": 5, "scientif": [5, 7], "comput": [2, 5, 7, 11], "base": [5, 7], "joint": 5, "work": 5, "kapidani": [5, 7], "codecasa": [5, 7], "book": 5, "design": 5, "provid": [5, 11], "introduct": 5, "exampl": 5, "galerkin": 5, "acoustiv": 5, "mix": [3, 5], "formul": [5, 6, 9, 10], "disconitinu": 5, "variant": 5, "approxim": [5, 10, 11], "conform": [5, 9], "each": [5, 8, 9, 10, 11], "other": [5, 10], "ansatz": 5, "featur": [5, 10], "differ": [5, 9, 10], "full": [5, 8], "mathemat": 5, "kcschoberl21": [5, 7], "wkcs23": [], "instal": 5, "maxwel": 5, "bernard": [5, 7], "lorenzo": [5, 7], "joachim": [5, 7], "sch": [5, 7], "\u00f6": [5, 7], "berl": [5, 7], "arbitrari": [5, 7], "depend": [5, 7, 9, 10, 11], "d": [3, 5, 6, 7, 8], "axwel": [5, 7], "phy": [5, 7], "433": [5, 7], "paper": [5, 7], "No": [5, 7], "110184": [5, 7], "2021": [5, 7], "doi": [5, 7], "1016": [5, 7], "jcp": [5, 7], "marku": [5, 7], "polynomi": [5, 7, 9, 11], "degre": [5, 7, 10], "acoust": [5, 7], "2023": [], "arxiv": 7, "2312": [], "14716": [], "section": [6, 8, 9, 10, 11], "describ": 6, "detail": 6, "how": 6, "appli": [6, 10, 11], "neglect": 6, "possibl": [6, 8, 10, 11], "here": [6, 11], "curl": [3, 6, 8], "text": [6, 10], "suitabl": [1, 6], "subset": 6, "mathbb": [6, 8, 9, 10], "close": 6, "constitut": 6, "relat": 6, "color": 6, "emph1": 6, "varepsilon": [6, 9, 10], "emph2": 6, "mu": [6, 9, 10], "permitt": 6, "permeabl": 6, "question": [6, 11], "weak": [6, 8, 9, 10], "eh": 6, "assum": 6, "homogen": 6, "abov": [6, 8, 9], "rewritten": 6, "377423e": [], "08": [1, 3], "270705e": [], "187766e": [], "version_info": [], "modulenotfounderror": [], "input": [], "name": 3, "9933350511855035": [], "1540": [], "06247129244141534": [], "381633e": [], "088": [], "160984e": [], "392659e": [], "35": [], "409019e": [], "322534e": [], "371604e": [], "0x7abf0dc794b0": [], "0x7abf0dc793f0": [], "0x7abf0dc78df0": [], "036663e": [], "3e": 2, "300": [], "dofss": 2, "when": [10, 11], "type": [], "problem": [8, 9, 11], "approach": [9, 10, 11], "typic": 10, "choic": [10, 11], "implicit": 10, "explicit": 11, "former": 10, "matrix": [8, 9, 10, 11, 12], "compos": [9, 10], "differenti": [1, 10, 12], "although": 10, "scheme": [10, 12], "uncondition": 10, "stabl": [3, 10, 12], "independ": 10, "factor": [10, 11], "larg": 10, "doe": 10, "well": [10, 12], "veri": [10, 11], "number": [9, 10], "freedom": 10, "hand": [8, 10], "mere": [8, 10], "There": 10, "exist": [9, 10], "sever": 10, "conveni": 10, "downsid": 10, "stabil": 10, "qualiti": 10, "finer": 10, "largest": 10, "guarante": 10, "popular": 10, "afdsj": [], "adf": [], "dg": 10, "appelogr20": [7, 10], "daniel": 7, "appel": 7, "fortino": 7, "garcia": 7, "olof": 7, "runborg": 7, "waveholtz": 7, "solut": [3, 7], "helmholtz": [7, 10], "siam": 7, "journal": [5, 7], "42": 7, "a1950": 7, "a1983": 7, "2020": [7, 10], "1137": 7, "19m1299062": 7, "apart": [10, 11], "obviou": [10, 11], "better": 10, "applic": [8, 10, 11], "solver": [7, 10], "interest": 10, "frequenc": 10, "scatter": [], "precondition": [7, 10], "stolk": 7, "eigenvalu": [7, 10, 11], "lite": 10, "pinvit": 10, "lobpcg": 10, "earli": 10, "2000": 10, "knyazev": 7, "filter": [7, 10], "2e5": 10, "hit": 10, "googl": 10, "scholar": 10, "2e4": 10, "sinc": [9, 10, 11], "2022": 10, "low": 10, "onli": [9, 10], "hexahedr": [9, 10], "grid": [9, 10], "5e4": 10, "6e3": 10, "numer": [10, 11], "flux": 10, "penalti": 10, "both": 10, "been": 10, "around": 10, "late": 10, "1960": 10, "1970": 10, "still": 10, "wide": 10, "engin": 10, "sto21": [7, 10], "christiaan": 7, "a3469": 7, "a3502": 7, "20m1359997": 7, "kny01": [7, 10], "andrew": 7, "toward": 7, "optim": [1, 7, 11], "precondit": 7, "eigensolv": 7, "local": [7, 9, 10], "conjug": 7, "gradient": 7, "517": 7, "541": 7, "2001": 7, "s1064827500366124": 7, "nw24": [7, 10], "lothar": 7, "nannen": 7, "krylov": 7, "2024": [5, 7], "2402": 7, "08515": 7, "idea": [9, 11], "stoke": 10, "theorem": 10, "my": [], "caption": [], "visual": [3, 9, 10], "lead": 10, "quantiti": 10, "interlac": 10, "pointwis": 10, "mimmick": 10, "l_": [9, 10], "k": [9, 10], "point": [1, 3, 8, 9, 10, 11], "valu": 10, "x_i": 10, "y_j": 10, "z_k": 10, "tangenti": [9, 10], "stagger": 10, "h_": 10, "tfrac": [8, 10], "approx": 10, "e_": 10, "satisfi": 10, "leap": [10, 12], "frog": [10, 12], "size": 10, "_": [9, 10], "frac": [8, 10, 12], "quad": [9, 10, 12], "carri": 10, "goal": 10, "gener": 10, "tetrahedr": [9, 10], "157738e": [], "basi": [1, 8, 10, 11], "focuss": 10, "go": 10, "dimens": [10, 11], "primal": [1, 9, 10, 11], "interpret": 10, "orang": 10, "circl": 10, "piecewis": 10, "constant": 10, "dark": 10, "grei": 10, "consist": [9, 10], "four": [9, 10], "proce": 10, "unknown": 10, "neighbour": [9, 10], "continu": [9, 10], "tild": [8, 10], "c_j": 10, "bigcup_": [9, 10], "c_": 10, "mathcal": [8, 10], "_j": 10, "w_": 10, "setr": [], "foral": [9, 10], "cont": 10, "v_": [10, 12], "global": [8, 9, 10], "w": 10, "zero": [9, 10, 11], "outsid": [9, 10], "origin": 10, "pose": 10, "sum": [], "rot": 10, "sum_": [8, 9, 10], "semi": [10, 12], "ultra": 10, "variat": 10, "remark": 10, "aris": 10, "across": 10, "part": 10, "help": 10, "skew": 10, "symmetri": 10, "immedi": 10, "conserv": 10, "simplex": 10, "motiv": 10, "construct": 11, "symbol": [1, 3], "correspond": [1, 11], "nodal": [1, 9, 11], "higher": 1, "sparsiti": [1, 11], "pattern": [1, 11], "328608e": 1, "596200e": 2, "352700e": 2, "824640e": 2, "0x7efeaab2eef0": 2, "0x7efeaab2f5f0": 2, "0x7efeaab2c9f0": 2, "879820000000851": 2, "030279e": 2, "sss": 2, "averag": [2, 3], "1546": [], "424357e": [], "home": [], "mwess": [], "lib": [], "python3": [], "site": [], "project": 8, "__init__": [], "py": 11, "63": [], "userwarn": [], "unabl": [], "axes3d": [], "multipl": [], "being": [], "pip": [], "3d": [0, 5], "warn": 3, "0x737388a1c9b0": [], "0x737388a1c270": [], "0x73738c7f48b0": [], "11327": [], "613625e": [], "309617e": [], "0x7c895f45ba70": [], "0x7c895d73d2b0": [], "0x7c895f343330": [], "552361e": [], "50": 1, "342434e": [], "498972e": [], "265443e": [], "0x7f00a4f15030": [], "0x7f00a4fa5530": [], "0x7f00a4fa5d30": [], "735248e": [], "274404e": [], "0x71581d5d7330": [], "0x71581d5d4330": [], "0x71581d5d5df0": [], "586595e": [], "408633e": [], "0x7aa7f4fbe670": [], "0x7aa7f4fbe0b0": [], "0x7aa7f4fbff70": [], "drawtim": 2, "drawnow": 2, "computetim": 2, "45": 9, "34500765800476": [], "553796e": [], "wkcs24": [5, 7], "physic": [5, 7, 8, 9], "page": [5, 7], "113196": [5, 7], "org": [5, 7], "decomposit": 9, "subdomain": 9, "bar": 9, "usual": [9, 11], "triangl": 9, "tetrahedra": 9, "therebi": 9, "barycent": 9, "connect": 9, "midpoint": 9, "edg": 9, "creat": 9, "quadrilater": 9, "share": 9, "vertex": 9, "triangular": 9, "split": 9, "becom": 9, "center": 9, "red": 9, "blue": 9, "uniqu": 9, "varphi_c": 9, "node": [9, 11], "x_0": [9, 11], "ldot": [9, 11], "unit": [9, 11], "interv": [9, 11], "tensor": [9, 11], "x_": 9, "xj": 9, "_p": [8, 9], "fulfil": [9, 12], "delta_": 9, "remap": 9, "onto": 9, "identifi": 9, "same": [8, 9, 11], "support": 9, "first": 9, "command": 9, "h1": [8, 9], "h1primalcel": [9, 11], "ne": 9, "take": [8, 9], "look": 9, "deform": [8, 9], "euler_angl": [3, 8, 9], "ones": [9, 11], "previou": 9, "e_k": 9, "e_0": 9, "e_1": 9, "canon": 9, "correct": 9, "covari": [8, 9], "transform": [8, 9, 11], "f_c": 9, "jacobian": [8, 9], "varphi": 9, "_c": 9, "definit": 9, "see": 9, "whole": [8, 9], "slightli": 9, "more": 9, "mesh_curl": 9, "435": 9, "hcurl": 9, "hcurldualcel": [8, 9, 11], "14": 9, "5714285714285716": 9, "gfu_curl": 9, "decompos": 9, "trace": [3, 9], "extens": 9, "minor": 9, "requir": [8, 9], "cube": [9, 11], "trilinear": 9, "semidiscret": 9, "written": 9, "down": 9, "similar": [8, 9, 11], "pleasant": 9, "check": 9, "By": 9, "face": [3, 9], "3p": 9, "12": 9, "6p": 9, "4p": 9, "verifi": [8, 9], "hold": 9, "unit_cub": [8, 9, 11], "hcurlprimalcel": [8, 9], "672": 9, "2736": 9, "7104": 9, "14640": 9, "count": 9, "nedg": 9, "nface": 9, "nelement": 9, "52": 9, "752": 9, "2964": 9, "7552": 9, "15380": 9, "would": 11, "chosen": 11, "classic": 11, "techniqu": 11, "x_p": 11, "yet": 11, "specifi": 11, "choos": 11, "exactli": 11, "quadratur": 11, "vanish": 11, "everi": 11, "except": 11, "expect": 11, "spars": 11, "gauss": 11, "radau": 11, "fix": 11, "treat": 11, "endpoint": 11, "gaussian": 11, "weight": 11, "easili": 11, "compar": 11, "exact": 11, "smooth": 11, "even": 11, "standard": 11, "fit": 11, "formula": 11, "trig_point": 11, "arrai": 11, "px": 11, "ob": 11, "spy": 11, "todens": 11, "seen": 11, "irs_f": 11, "m_diag": 11, "deletezeroel": [3, 11], "coupl": 11, "entri": 11, "fes_curl": 11, "irs_fes_curl": 11, "m_lump": 11, "store": [8, 11], "togeth": 11, "structur": 11, "less": 11, "exploit": 11, "access": 11, "m_exact": 11, "exacttim": 11, "m_superspars": 11, "stime": 11, "superspars": 11, "m_exact_inv": 11, "sparsecholeski": [3, 11], "m_supersparse_inv": 11, "tmp2": 11, "trialfunct": 8, "dx_vol": 8, "dx_edg": 8, "cross": [3, 8], "geometr": 8, "contribut": 8, "cancel": 8, "out": 8, "permut": 8, "equival": 8, "class": 8, "reduct": 8, "memori": 8, "cost": 8, "realiz": 8, "flag": 8, "curl_gf": 8, "curlt": 8, "gft": 8, "geometry_fre": 8, "100": [3, 8], "miss": 8, "spacial": 8, "nameerror": [], "befor": 8, "launch": 8, "hdiv": 8, "hdivdualcel": [], "h1primal": 8, "hdivdual": 8, "hiv": [], "side": 8, "mass_h1_inv": 8, "dx_h1": 8, "rh": 8, "gfp": 8, "1552574634552002": [], "014244794845581055": [], "16961169242858887": [], "029892921447753906": [], "15879392623901367": [], "013952970504760742": [], "15647482872009277": [], "024718046188354492": [], "attributeerror": [], "attribut": [], "15273332595825195": [], "014251708984375": [], "22423505783081055": [], "03263998031616211": [], "hprimaldualcel": [], "1410810947418213": [], "014180183410644531": [], "1572096347808838": [], "02831745147705078": [], "13854": [], "15579557418823242": [], "014397859573364258": [], "18164634704589844": [], "031702280044555664": [], "30510": [], "runtimeerror": [], "file": [], "401": [], "obj": [], "show": 3, "arg": [], "kwarg": [], "399": [], "402": [], "kwargs_with_default": [], "width": [], "height": [], "403": [], "404": [], "filenam": [], "405": [], "generatehtml": [], "widget": [], "76": [], "basewebguiscen": 3, "self": [], "74": [], "webguiwidget": [], "layout": [], "75": [], "encod": [], "getdata": [], "77": [], "displai": [], "260": [], "webglscen": [], "set_minmax": [], "257": [], "encodedata": [], "dtype": [], "float32": [], "259": [], "none": [], "rais": [], "cannot": [], "typ": [], "262": [], "263": [], "la": [], "basevector": [], "1414651870727539": [], "013903617858886719": [], "18274283409118652": [], "025954246520996094": [], "14151358604431152": [], "013034343719482422": [], "15631484985351562": [], "02606368064880371": [], "154870": 8, "344250": 8, "14334511756896973": [], "014346122741699219": [], "1761927604675293": [], "026534080505371094": [], "40": 8, "16558480262756348": [], "014064788818359375": [], "25002050399780273": [], "02897024154663086": [], "15956687927246094": [], "014201641082763672": [], "2367711067199707": [], "027669191360473633": [], "90": [], "13974499702453613": [], "01331949234008789": [], "1559600830078125": [], "024106740951538086": [], "150": 8, "14138245582580566": [], "013830423355102539": [], "20544767379760742": [], "028992652893066406": [], "1550123691558838": [], "014713525772094727": [], "22198820114135742": [], "026645660400390625": [], "includ": 8, "distribut": 8, "1411738395690918": [], "014023542404174805": [], "1531846523284912": [], "02953052520751953": [], "_h": 8, "int_t": 8, "mass_hdiv_inv": 8, "15952825546264648": [], "016004562377929688": [], "21116304397583008": [], "036969661712646484": [], "14396190643310547": [], "014764785766601562": [], "15488696098327637": [], "024342060089111328": [], "14777612686157227": [], "014744281768798828": [], "24978208541870117": [], "034796953201293945": [], "15024352073669434": [], "014367341995239258": [], "214324951171875": [], "03016209602355957": [], "m_u": [], "dot": [], "m_v": 12, "v_0": 12, "u_0": [], "u_": [], "u_j": [], "14087343215942383": [], "013223409652709961": [], "1841111183166504": [], "03317904472351074": [], "1401073932647705": [], "013926029205322266": [], "15143561363220215": [], "024031400680541992": [], "m_p": 12, "p_0": 12, "p_": 12, "p_j": 12, "13863372802734375": [], "014627456665039062": [], "2688016891479492": [], "02542257308959961": [], "known": 12, "sigma": 12, "spectral": 12, "radiu": [3, 12], "1457195281982422": [], "013231277465820312": [], "15278911590576172": [], "04138064384460449": [], "parital_t": [], "1678910255432129": [], "015158414840698242": [], "156968355178833": [], "029218196868896484": [], "14244389533996582": [], "013953447341918945": [], "2671644687652588": [], "03650474548339844": [], "14476919174194336": [], "013946533203125": [], "15114641189575195": [], "027591228485107422": [], "14116525650024414": [], "014034509658813477": [], "151580810546875": [], "02973484992980957": [], "locnr": [], "105662": [], "394338": [], "894338": [], "605662": [], "803561": [], "0982194": [], "0536948": [], "555556": [], "0778847": [], "0416667": [], "36656": [], "418661": [], "290669": [], "0296385": [], "725312": [], "0915628": [], "00999851": [], "51153": [], "0733767": [], "00616387": [], "341717": [], "391248": [], "0610607": [], "273846": [], "00365801": [], "316355": [], "227882": [], "0022025": [], "65232": [], "1840486526489258": [], "10226964950561523": [], "141755104064941": [], "0056302547454833984": [], "18722081184387207": [], "0008227825164794922": [], "14056038856506348": [], "014088630676269531": [], "15216398239135742": [], "030299901962280273": [], "1457281112670898": [], "0982968807220459": [], "94332766532898": [], "0051441192626953125": [], "1732780933380127": [], "000850677490234375": [], "341424e": [], "0x72bfd6786630": [], "0x72bfd4b248b0": [], "0x72bfd68f50b0": [], "44": [], "62706017494202": [], "707495e": [], "14597392082214355": [], "01627063751220703": [], "15703964233398438": [], "025795698165893555": [], "2047088146209717": [], "11965441703796387": [], "138391017913818": [], "005194664001464844": [], "17090559005737305": [], "0008075237274169922": [], "16479754447937012": [], "014780998229980469": [], "15702223777770996": [], "02642369270324707": [], "1339991092681885": [], "10081267356872559": [], "03842544555664": [], "0051000118255615234": [], "18155431747436523": [], "0008234977722167969": [], "1769404411315918": [], "020380258560180664": [], "17107176780700684": [], "024678468704223633": [], "122352": [], "913511037826538": [], "17615032196044922": [], "74958562850952": [], "008779048919677734": [], "4080026149749756": [], "003787517547607422": [], "14220356941223145": [], "013993024826049805": [], "2083144187927246": [], "03446626663208008": [], "0105292797088623": [], "21221518516540527": [], "43076729774475": [], "007737159729003906": [], "389972448348999": [], "0016677379608154297": [], "384611e": [], "0x7af0b864c4f0": [], "0x7af10877e570": [], "0x7af10877ecb0": [], "11096501350403": [], "603363e": [], "14337658882141113": [], "016633987426757812": [], "15613245964050293": [], "02317190170288086": [], "3169796466827393": [], "24944090843200684": [], "955610036849976": [], "007941484451293945": [], "39960145950317383": [], "0014188289642333984": [], "382613e": [], "0x7e7b5c13bcf0": [], "0x7e7b5dd9fbf0": [], "0x7e7b5dd9f0b0": [], "98267078399658": [], "421309e": [], "1393754482269287": [], "015882492065429688": [], "1510932445526123": [], "02279043197631836": [], "2654900550842285": [], "2561802864074707": [], "194032430648804": [], "008639097213745117": [], "3899693489074707": [], "0015819072723388672": [], "dualmesh": [], "gui": 3, "lx": 3, "ly": 3, "lz": 3, "obstacl": 3, "rad": 3, "experiment": 3, "t0": 3, "occ": 3, "wg": 3, "box": 3, "pnt": 3, "z": 3, "inflow": 3, "sphere": 3, "occgeometri": 3, "clip": 3, "120": 3, "15": 3, "fes_": 3, "hcurldualcells3d": 3, "fes_h": 3, "hcurlprimalcells3d": 3, "gf": 3, "dxe": 3, "dse": 3, "bf_mix": 3, "1280304": 3, "prepar": 3, "massh_inv": 3, "bfm_e": 3, "masse_inv": 3, "masse_surf": 3, "masse_surf_inv": 3, "getdof": 3, "inconsist": 3, "silenc": 3, "check_unus": 3, "plx": 3, "lfr": 3, "rhsefunc": 3, "lambda": 3, "01": 3, "scenebd": [], "draw_surf": 3, "sceneh": 3, "altshap": 3, "17777777777777762": [], "831430e": [], "35555555555555435": [], "598219e": [], "5333333333333311": [], "887032e": [], "7111111111111078": [], "012115e": [], "8888888888888845": [], "500": [], "951472e": [], "0666666666666653": [], "600": [], "115749e": [], "2444444444444531": [], "700": [], "067154e": [], "422222222222241": [], "800": [], "005547e": [], "6000000000000287": [], "900": [], "294784e": [], "7777777777778165": [], "093375e": [], "9555555555556043": [], "1100": [], "110797e": [], "133333333333392": [], "1200": [], "211048e": [], "3111111111111797": [], "1300": [], "056489e": [], "4888888888889675": [], "1400": [], "205722e": [], "6666666666667553": [], "1500": [], "146143e": [], "844444444444543": [], "1600": [], "149856e": [], "022222222222331": [], "1700": [], "331378e": [], "2000000000001187": [], "1800": [], "118499e": [], "3777777777779066": [], "1900": [], "036926e": [], "5555555555556944": [], "357804e": [], "733333333333482": [], "2100": [], "093240e": [], "91111111111127": [], "2200": [], "196624e": [], "19479799270629883": [], "02530384063720703": [], "2621617317199707": [], "06622910499572754": [], "gfe_anim": 3, "255827e": [], "ss": [], "13": [], "18": [], "19": [], "finsh": 3, "14386940002441406": [], "01591348648071289": [], "15303897857666016": [], "02896285057067871": [], "2036757469177246": [], "027491092681884766": [], "23320865631103516": [], "07455897331237793": [], "081141e": [], "18926262855529785": [], "02477884292602539": [], "25863194465637207": [], "08016109466552734": [], "eval": 3, "042641e": [], "1747128963470459": [], "023659229278564453": [], "20588016510009766": [], "05399680137634277": [], "749302e": 3, "1346902847290039": [], "013793468475341797": [], "17264318466186523": [], "023189544677734375": [], "13655972480773926": 8, "014809608459472656": 8, "16208624839782715": 8, "028594493865966797": 8}, "objects": {}, "objtypes": {}, "objnames": {}, "titleterms": {"exampl": 0, "A": 1, "plane": 1, "wave": 1, "squar": 1, "ring": 2, "reson": [2, 10], "gener": [2, 9], "geometri": 2, "pml": 2, "sourc": 2, "coeffici": 2, "space": [2, 9, 10], "bilinear": 2, "form": 2, "time": [2, 6, 12], "loop": 2, "instal": 4, "get": 4, "ngsolv": [4, 5], "dual": [4, 5, 6, 9, 10], "cell": [4, 5, 6, 9, 10], "method": [4, 5, 6, 10], "add": 4, "us": 4, "pip": 4, "work": 4, "code": 4, "test": 4, "dualcellspac": 4, "The": [5, 6, 8, 9], "tabl": 5, "content": 5, "refer": 5, "domain": 6, "maxwel": 6, "system": 6, "problem": [6, 10], "set": [6, 10], "basic": 10, "idea": 10, "state": 10, "art": 10, "from": 10, "fdtd": 10, "construct": [9, 10], "discret": [8, 9], "mass": 11, "lump": 11, "step": 12, "introduct": 10, "reason": 10, "explicit": 10, "scatter": 10, "type": 10, "galerkin": 10, "high": 10, "order": 10, "barycentr": 9, "mesh": 9, "refin": 9, "2d": 9, "map": 9, "basi": 9, "function": 9, "x_p": 9, "mathrm": 9, "grad": 9, "mathcal": 9, "t": 9, "tild": 9, "curl": 9, "summeri": 9, "three": 9, "dimens": 9, "assembl": [], "differenti": 8, "oper": 8, "etud": [], "comput": 8, "gradient": 8, "gaussian": 8, "peak": 8, "\u00e9tude": 8, "dualmesh": [], "an": [], "electromagnet": 3, "waveguid": 3, "3d": 3, "bibliographi": 7}, "envversion": {"sphinx.domains.c": 3, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 9, "sphinx.domains.index": 1, "sphinx.domains.javascript": 3, "sphinx.domains.math": 2, "sphinx.domains.python": 4, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx.ext.intersphinx": 1, "sphinxcontrib.bibtex": 9, "sphinx": 60}, "alltitles": {"Examples": [[0, "examples"]], "A plane wave on a square": [[1, "a-plane-wave-on-a-square"]], "Ring resonator": [[2, "ring-resonator"]], "Generating the geometry": [[2, "generating-the-geometry"]], "PML geometry": [[2, "pml-geometry"]], "Source": [[2, "source"]], "Coefficients": [[2, "coefficients"]], "Spaces": [[2, "spaces"]], "Bilinear forms": [[2, "bilinear-forms"]], "Time loop": [[2, "time-loop"]], "Electromagnetic 3d waveguide": [[3, "electromagnetic-3d-waveguide"]], "Installation": [[4, "installation"]], "Get NGSolve": [[4, "get-ngsolve"]], "Install the Dual Cell Method add-on": [[4, "install-the-dual-cell-method-add-on"]], "Using pip": [[4, "using-pip"]], "Working with the code": [[4, "working-with-the-code"]], "Test the dualcellspaces installation": [[4, "test-the-dualcellspaces-installation"]], "The Dual Cell Method in NGSolve": [[5, "the-dual-cell-method-in-ngsolve"]], "Table of Contents": [[5, "table-of-contents"]], "References": [[5, "references"]], "The dual cell method for the time-domain Maxwell system": [[6, "the-dual-cell-method-for-the-time-domain-maxwell-system"]], "Problem setting": [[6, "problem-setting"]], "Bibliography": [[7, "bibliography"]], "The discrete differential operators": [[8, "the-discrete-differential-operators"]], "\u00c9tude: computing the gradient of a Gaussian peak": [[8, "etude-computing-the-gradient-of-a-gaussian-peak"]], "Discrete spaces on dual cells": [[9, "discrete-spaces-on-dual-cells"]], "Barycentric mesh refinement in 2d": [[9, "barycentric-mesh-refinement-in-2d"]], "Mapped basis functions": [[9, "mapped-basis-functions"]], "The space X_P^{\\mathrm{grad}}(\\mathcal T)": [[9, "the-space-x-p-mathrm-grad-mathcal-t"]], "The space \\tilde X_P^{\\mathrm{curl}}(\\tilde{\\mathcal T})": [[9, "the-space-tilde-x-p-mathrm-curl-tilde-mathcal-t"]], "Summery of the construction": [[9, "summery-of-the-construction"]], "Generalization to three dimensions": [[9, "generalization-to-three-dimensions"]], "Introduction": [[10, "introduction"]], "Reasons for explicit methods": [[10, "reasons-for-explicit-methods"]], "scattering type problems": [[10, "scattering-type-problems"]], "resonance type problems": [[10, "resonance-type-problems"]], "State of the art": [[10, "state-of-the-art"]], "From FDTD to the dual cell method": [[10, "from-fdtd-to-the-dual-cell-method"]], "Basic idea of the dual cell construction": [[10, "basic-idea-of-the-dual-cell-construction"]], "Galerkin setting and high order spaces": [[10, "galerkin-setting-and-high-order-spaces"]], "Mass lumping": [[11, "mass-lumping"]], "Time stepping": [[12, "time-stepping"]]}, "indexentries": {}})