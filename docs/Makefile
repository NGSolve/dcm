DEPENDS=intro.md _config.yml _toc.yml installation.md createimages.pdf maxwell.md maxwell/basic_ideas.md#first_steps.md \
.PHONY: clean

_build: $(DEPENDS)
	jupyter-book build --all .

createimages.pdf: createimages.tex 
	lualatex --shell-escape $< 
