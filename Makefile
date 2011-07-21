PYTHON = python
COVERAGE_REPORT=coverage

clean:
	$(MAKE) -C doc clean
	-rm -rf build
	-rm *-stamp

distclean: clean
	-rm MANIFEST
	-rm $(COVERAGE_REPORT)
	@find . -name '*.py[co]' \
		 -o -name '*.a' \
		 -o -name '*,cover' \
		 -o -name '.coverage' \
		 -o -iname '*~' \
		 -o -iname '*.kcache' \
		 -o -iname '*.pstats' \
		 -o -iname '*.prof' \
		 -o -iname '#*#' | xargs -L10 rm -f

# Print out info for possible install methods
check-version-info:
	$(PYTHON) -c 'from nisext.testers import info_from_here; info_from_here("formula")'

# Run tests from installed code
installed-tests:
	$(PYTHON) -c 'from nisext.testers import tests_installed; tests_installed("formula")'

# Run tests from installed code
sdist-tests:
	$(PYTHON) -c 'from nisext.testers import sdist_tests; sdist_tests("formula")'

bdist-egg-tests:
	$(PYTHON) -c 'from nisext.testers import bdist_egg_tests; bdist_egg_tests("formula")'

source-release: distclean
	python -m compileall .
	make distclean
	python setup.py sdist --formats=gztar,zip
