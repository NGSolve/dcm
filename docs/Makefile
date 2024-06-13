DEPENDS=intro.md _config.yml _toc.yml installation.md createimages.pdf #first_steps.md \
.PHONY: clean

_build: $(DEPENDS)
	jupyter-book build --all .

createimages.pdf: createimages.tex 
	lualatex --shell-escape $< 
