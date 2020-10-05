.PHONY: info build shell py2 py3 test ftest publish

info:
	@echo ""
	@echo "Available Make Commands"
	@echo ""
	@echo "build: builds the image"
	@echo "shell: starts a bash shell in the container
	@echo "py2:	maps ~/Documents/nd2s to /var/nd2s and runs a Python 2.7 interpreter"
	@echo "py3:	maps ~/Documents/nd2s to /var/nd2s and runs a Python 3.4 interpreter"
	@echo "test: runs all unit tests (in Python 3.4)"
	@echo "ftest: runs all functional tests (requires specific ND2 files that are not publicly available"
	@echo "publish: publishes the code base to PyPI (maintainers only)"
	@echo ""

build:
	docker build -t jimrybarski/nd2reader .

shell:
	xhost local:root; docker run --rm -v $(CURDIR):/opt/nd2reader -v ~/nd2s:/var/nd2s -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=unix$(DISPLAY) -it jimrybarski/nd2reader bash

py2:
	xhost local:root; docker run --rm -v $(CURDIR):/opt/nd2reader -v ~/nd2s:/var/nd2s -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=unix$(DISPLAY) -it jimrybarski/nd2reader python2.7

py3:
	xhost local:root; docker run --rm -v $(CURDIR):/opt/nd2reader -v ~/nd2s:/var/nd2s -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=unix$(DISPLAY) -it jimrybarski/nd2reader python3.4

test:	build
	docker run --rm -v $(CURDIR):/opt/nd2reader -it jimrybarski/nd2reader python3.4 test.py
	docker run --rm -v $(CURDIR):/opt/nd2reader -it jimrybarski/nd2reader python2.7 test.py

ftest:	build
	xhost local:root; docker run --rm -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=unix$(DISPLAY) -v $(CURDIR):/opt/nd2reader -v ~/nd2s:/var/nd2s -it jimrybarski/nd2reader python3.4 /opt/nd2reader/ftest.py
	xhost local:root; docker run --rm -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=unix$(DISPLAY) -v $(CURDIR):/opt/nd2reader -v ~/nd2s:/var/nd2s -it jimrybarski/nd2reader python2.7 /opt/nd2reader/ftest.py
	
publish:
	python setup.py sdist upload -r pypi
