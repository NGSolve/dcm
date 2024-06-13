DEPENDS=intro.md first_steps.md _config.yml _toc.yml installation.md
.PHONY: clean

_build: $(DEPENDS)
	jupyter-book build --all .
