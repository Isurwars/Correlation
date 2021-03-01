SUBDIRS = src
dist_doc_DATA = README.md LICENSE

TST_DIR = test

EXE = $(top_builddir)/bin/correlation$(EXEEXT)

tests: test1 test2 test3 test4 test5

test1: $(TST_DIR)/test_1/xSi.car
	$(EXE) $(TST_DIR)/test_1/xSi.car

test2: $(TST_DIR)/test_2/Graphene.car
	$(EXE) $(TST_DIR)/test_2/Graphene.car

test3: $(TST_DIR)/test_3/aPd.cell
	$(EXE) $(TST_DIR)/test_3/aPd.cell

test4: $(TST_DIR)/test_4/aPdH.dat
	$(EXE) -i $(TST_DIR)/test_4/aPdH.txt $(TST_DIR)/test_4/aPdH.dat

test5: $(TST_DIR)/test_5/lBi.car
	$(EXE) $(TST_DIR)/test_5/lBi.car

test6: $(TST_DIR)/test_6/aPdH_3144.car
	$(EXE) $(TST_DIR)/test_6/aPdH_3144.car


clean-local:
	@rm config.status config.log configure~
	@rm -r autom4te.cache/
	@rm aclocal.m4 Makefile
	@rm src/Makefile
	@rm src/*.o
	@rm -r src/.deps/
	@rm bin/.dirstamp
