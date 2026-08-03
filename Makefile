# Makefile for ProPrep - AmberTools integration
#
# When built as part of AmberTools, config.h provides:
#   $(PYTHON)  - path to the Python interpreter
#   $(BINDIR)  - path to $AMBERHOME/bin
#
# For standalone use:
#   make install PYTHON=python3 AMBERHOME=/path/to/amber

-include $(AMBERHOME)/AmberTools/src/config.h

# Fallback if not building within AmberTools
PYTHON ?= python3
BINDIR ?= $(AMBERHOME)/bin

install:
	$(PYTHON) setup.py install --prefix=$(AMBERHOME) --install-scripts=$(BINDIR)

# Release guard: verify the four version pins are in lockstep before building.
# Optionally pass VERSION=X.Y.Z to also assert the agreed value.
check-version:
	bash scripts/check_version_lockstep.sh $(VERSION)

clean:
	-rm -rf build/ dist/ *.egg-info src/*.egg-info
	-find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
	-find . -name "*.pyc" -delete 2>/dev/null || true

uninstall:
	@echo "To uninstall, remove proprep from:"
	@echo "  $(AMBERHOME)/lib/python*/site-packages/proprep/"
	@echo "  $(BINDIR)/proprep"
	@echo "  $(BINDIR)/mpsa"
