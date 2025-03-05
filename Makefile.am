SUBDIRS = src
dist_doc_DATA = README.md LICENSE

TST_DIR = $(top_builddir)/test
BIN_DIR = $(top_builddir)/bin

EXE = $(BIN_DIR)/correlation$(EXEEXT)


###################################################################
## Test options with automatic testing of results
###################################################################
tests: test1 test2 test3 test4 test5 test6

test1:
	@echo "--------------------------------------------------"
	@echo " Running test 1: (Crystalline Silicon)"
	@echo "--------------------------------------------------"
	$(EXE) -a 0.1 $(TST_DIR)/test_1/xSi.car
	@diff --strip-trailing-cr --strip-trailing-cr $(TST_DIR)/test_1/xSi_g.csv $(TST_DIR)/test_1/g.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo "----------------------------------------------------------"; \
		echo " Test 1: PDF (g) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo "  Warning: PDF (g) different from expected."; \
		echo "--------------------------------------------------"; \
	fi
	@diff --strip-trailing-cr --strip-trailing-cr $(TST_DIR)/test_1/xSi_J.csv $(TST_DIR)/test_1/J.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo " Test 1: RDF (J) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo " Warning: RDF (J) different from expected."; \
		echo "----------------------------------------------------------"; \
	fi
	@diff --strip-trailing-cr --strip-trailing-cr $(TST_DIR)/test_1/xSi_G_.csv $(TST_DIR)/test_1/G_.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo " Test 1: rPDF (G) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo " Warning: rPDF (G) different from expected."; \
		echo "----------------------------------------------------------"; \
	fi
	@diff --strip-trailing-cr --strip-trailing-cr $(TST_DIR)/test_1/xSi_Z.csv $(TST_DIR)/test_1/Z.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo " Test 1: Coordination Number (Z) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo " Warning: Coordination Number (Z) different from expected."; \
		echo "----------------------------------------------------------"; \
	fi
	@diff --strip-trailing-cr --strip-trailing-cr $(TST_DIR)/test_1/xSi_PAD.csv $(TST_DIR)/test_1/PAD.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo " Test 1: PAD (F) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo " Warning: PAD (F) different from expected."; \
		echo "----------------------------------------------------------"; \
	fi

test2:
	@echo "--------------------------------------------------"
	@echo " Running test 2 : (Graphene Single-Layer)"
	@echo "--------------------------------------------------"
	$(EXE) -b 1.3 $(TST_DIR)/test_2/Graphene.car
	@diff --strip-trailing-cr $(TST_DIR)/test_2/Graphene_g.csv $(TST_DIR)/test_2/g.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo "----------------------------------------------------------"; \
		echo " Test 2: PDF (g) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo "  Warning: PDF (g) different from expected."; \
		echo "--------------------------------------------------"; \
	fi
	@diff --strip-trailing-cr $(TST_DIR)/test_2/Graphene_J.csv $(TST_DIR)/test_2/J.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo " Test 2: RDF (J) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo " Warning: RDF (J) different from expected."; \
		echo "----------------------------------------------------------"; \
	fi
	@diff --strip-trailing-cr --strip-trailing-cr $(TST_DIR)/test_2/Graphene_G_.csv $(TST_DIR)/test_2/G_.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo " Test 2: rPDF (G) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo " Warning: rPDF (G) different from expected."; \
		echo "----------------------------------------------------------"; \
	fi
	@diff --strip-trailing-cr $(TST_DIR)/test_2/Graphene_Z.csv $(TST_DIR)/test_2/Z.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo " Test 2: Coordination Number (Z) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo " Warning: Coordination Number (Z) different from expected."; \
		echo "----------------------------------------------------------"; \
	fi
	@diff --strip-trailing-cr $(TST_DIR)/test_2/Graphene_PAD.csv $(TST_DIR)/test_2/PAD.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo " Test 2: PAD (F) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo " Warning: PAD (F) ifferent from expected."; \
		echo "----------------------------------------------------------"; \
	fi


test3:
	@echo "--------------------------------------------------"
	@echo " Running test 3 : (Amorphous Palladium)"
	@echo "--------------------------------------------------"
	$(EXE) -q 0.1 -S $(TST_DIR)/test_3/aPd.cell
	@diff --strip-trailing-cr $(TST_DIR)/test_3/aPd_g.csv $(TST_DIR)/test_3/g.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo "----------------------------------------------------------"; \
		echo " Test 3: PDF (g) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo "  Warning: PDF (g) different from expected."; \
		echo "--------------------------------------------------"; \
	fi
	@diff --strip-trailing-cr $(TST_DIR)/test_3/aPd_J.csv $(TST_DIR)/test_3/J.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo " Test 3: RDF (J) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo " Warning: RDF (J) different from expected."; \
		echo "----------------------------------------------------------"; \
	fi
	@diff --strip-trailing-cr $(TST_DIR)/test_3/aPd_G_.csv $(TST_DIR)/test_3/G_.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo " Test 3: rPDF (G) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo " Warning: rPDF (G) different from expected."; \
		echo "----------------------------------------------------------"; \
	fi
	@diff --strip-trailing-cr $(TST_DIR)/test_3/aPd_Z.csv $(TST_DIR)/test_3/Z.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo " Test 3: Coordination Number (Z) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo " Warning: Coordination Number (Z) different from expected."; \
		echo "----------------------------------------------------------"; \
	fi
	@diff --strip-trailing-cr $(TST_DIR)/test_3/aPd_PAD.csv $(TST_DIR)/test_3/PAD.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo " Test 3: PAD (F) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo " Warning: PAD (F) different from expected."; \
		echo "----------------------------------------------------------"; \
	fi


test4:
	@echo "--------------------------------------------------"
	@echo " Running test 4 : (Amorphous Palladium Hydride)"
	@echo "--------------------------------------------------"
	$(EXE) -i $(TST_DIR)/test_4/aPdH.txt $(TST_DIR)/test_4/aPdH.dat
	@diff --strip-trailing-cr $(TST_DIR)/test_4/aPdH_g.csv $(TST_DIR)/test_4/g.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo "----------------------------------------------------------"; \
		echo " Test 4: PDF (g) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo "  Warning: PDF (g) different from expected."; \
		echo "--------------------------------------------------"; \
	fi
	@diff --strip-trailing-cr $(TST_DIR)/test_4/aPdH_J.csv $(TST_DIR)/test_4/J.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo " Test 4: RDF (J) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo " Warning: RDF (J) different from expected."; \
		echo "----------------------------------------------------------"; \
	fi
	@diff --strip-trailing-cr $(TST_DIR)/test_4/aPdH_G_.csv $(TST_DIR)/test_4/G_.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo " Test 4: rPDF (G) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo " Warning: rPDF (G) different from expected."; \
		echo "----------------------------------------------------------"; \
	fi
	@diff --strip-trailing-cr $(TST_DIR)/test_4/aPdH_Z.csv $(TST_DIR)/test_4/Z.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo " Test 4: Coordination Number (Z) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo " Warning: Coordination Number (Z) different from expected."; \
		echo "----------------------------------------------------------"; \
	fi
	@diff --strip-trailing-cr $(TST_DIR)/test_4/aPdH_PAD.csv $(TST_DIR)/test_4/PAD.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo " Test 4: PAD (F) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo " Warning: PAD (F) different from expected."; \
		echo "----------------------------------------------------------"; \
	fi


test5:
	@echo "--------------------------------------------------"
	@echo " Running test 5 : (Liquid Bismuth)\n"
	@echo "--------------------------------------------------"
	$(EXE) -q 0.02 $(TST_DIR)/test_5/lBi.car
	@diff --strip-trailing-cr $(TST_DIR)/test_5/lBi_g.csv $(TST_DIR)/test_5/g.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo "----------------------------------------------------------"; \
		echo " Test 5: PDF (g) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo "  Warning: PDF (g) different from expected."; \
		echo "--------------------------------------------------"; \
	fi
	@diff --strip-trailing-cr $(TST_DIR)/test_5/lBi_J.csv $(TST_DIR)/test_5/J.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo " Test 5: RDF (J) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo " Warning: RDF (J) different from expected."; \
		echo "----------------------------------------------------------"; \
	fi
	@diff --strip-trailing-cr $(TST_DIR)/test_5/lBi_G_.csv $(TST_DIR)/test_5/G_.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo " Test 5: rPDF (G) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo " Warning: rPDF (G) different from expected."; \
		echo "----------------------------------------------------------"; \
	fi
	@diff --strip-trailing-cr $(TST_DIR)/test_5/lBi_Z.csv $(TST_DIR)/test_5/Z.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo " Test 5: Coordination Number (Z) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo " Warning: Coordination Number (Z) different from expected."; \
		echo "----------------------------------------------------------"; \
	fi
	@diff --strip-trailing-cr $(TST_DIR)/test_5/lBi_PAD.csv $(TST_DIR)/test_5/PAD.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo " Test 5: PAD (F) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo " Warning: PAD (F) different from expected."; \
		echo "----------------------------------------------------------"; \
	fi


test6:
	@echo "--------------------------------------------------"
	@echo " Running test 6 : (3000+ atoms)\n"
	@echo "--------------------------------------------------"
	$(EXE) -R 30.0 $(TST_DIR)/test_6/aPdH_3144.car
	@diff --strip-trailing-cr $(TST_DIR)/test_6/aPdH_3144_g.csv $(TST_DIR)/test_6/g.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo "----------------------------------------------------------"; \
		echo " Test 6: PDF (g) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo "  Warning: PDF (g) different from expected."; \
		echo "--------------------------------------------------"; \
	fi
	@diff --strip-trailing-cr $(TST_DIR)/test_6/aPdH_3144_J.csv $(TST_DIR)/test_6/J.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo " Test 6: RDF (J) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo " Warning: RDF (J) different from expected."; \
		echo "----------------------------------------------------------"; \
	fi
	@diff --strip-trailing-cr $(TST_DIR)/test_6/aPdH_3144_G_.csv $(TST_DIR)/test_6/G_.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo " Test 6: rPDF (G) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo " Warning: rPDF (G) different from expected."; \
		echo "----------------------------------------------------------"; \
	fi
	@diff --strip-trailing-cr $(TST_DIR)/test_6/aPdH_3144_Z.csv $(TST_DIR)/test_6/Z.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo " Test 6: Coordination Number (Z) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo " Warning: Coordination Number (Z) different from expected."; \
		echo "----------------------------------------------------------"; \
	fi
	@diff --strip-trailing-cr $(TST_DIR)/test_6/aPdH_3144_PAD.csv $(TST_DIR)/test_6/PAD.csv; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo " Test 6: PAD (F) successfully calculated."; \
		echo "----------------------------------------------------------"; \
	else \
		echo "----------------------------------------------------------"; \
		echo " Warning: PAD (F) different from expected."; \
		echo "----------------------------------------------------------"; \
	fi


test7:
	@echo "--------------------------------------------------"
	@echo " Stress Test : (Amorphous Palladium Deuterium)"
	@echo "--------------------------------------------------"
	$(EXE) $(TST_DIR)/test_7/aPdD-stress.car
	@echo "--------------------------------------------------"
	@echo " Stress Test successfully completed."
	@echo "--------------------------------------------------"

test8:
	@echo "--------------------------------------------------"
	@echo " Stress Test 2: (Amorphous Palladium Gold Hydrogen)"
	@echo "--------------------------------------------------"
	$(EXE) -r 40.0 $(TST_DIR)/test_8/aPdAuH-stress.car
	@echo "--------------------------------------------------"
	@echo " Stress Test 2 successfully completed."
	@echo "--------------------------------------------------"


###################################################################
## Clean testing mess to match exactly distributed files
###################################################################
clean-test:
	@rm -f $(TST_DIR)/test_1/xSi*.csv
	@rm -f $(TST_DIR)/test_2/Graphene*.csv
	@rm -f $(TST_DIR)/test_3/aPd*.csv
	@rm -f $(TST_DIR)/test_4/aPdH*.csv
	@rm -f $(TST_DIR)/test_5/lBi*.csv
	@rm -f $(TST_DIR)/test_6/aPdH_3144*.csv
	@rm -f $(TST_DIR)/test_7/aPdD-stress*.csv


###################################################################
## Clean compilation mess while keeping working testing files
###################################################################
clean-local: clean-test
	@rm -f config.status config.log configure~
	@rm -rf autom4te.cache/
	@rm -f aclocal.m4
	@rm -f src/Makefile
	@rm -f src/*.o
	@rm -rf src/.deps/
	@rm -f bin/.dirstamp

###################################################################
## Create pristine directories to match exactly distributed files
###################################################################
clean-all: distclean
	@rm -f config.status config.log configure~
	@rm -rf autom4te.cache/
	@rm -f aclocal.m4 Makefile
	@rm -f src/Makefile
	@rm -f src/*.o
	@rm -rf src/.deps/
	@rm -f bin/.dirstamp
