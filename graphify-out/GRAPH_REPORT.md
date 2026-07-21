# Graph Report - Correlation  (2026-07-20)

## Corpus Check
- 279 files · ~5,125,080 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 3937 nodes · 7175 edges · 256 communities (213 shown, 43 thin omitted)
- Extraction: 92% EXTRACTED · 8% INFERRED · 0% AMBIGUOUS · INFERRED: 558 edges (avg confidence: 0.8)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `6f7ca20d`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- PDF Generation Utilities
- CLI Option Validation Tests
- ChiralityCalculator.cpp
- PresetManager.cpp
- TrajectoryAnalyzer
- Atom
- KernelGenerationParams
- XYZReader
- LinearAlgebra.hpp
- vector
- AppBackend
- SIMDUtils.hpp
- Cellulose Example Data
- DistributionFunctions
- DistributionFunctions.cpp
- ThreadLocalDistances
- CliOptions
- PlotController
- InputValidator.cpp
- TEST_F
- AppBackend.cpp
- TEST_F
- CastepMdReader
- AnalysisSettings
- AppDefaults
- CifReader.cpp
- TEST_F
- FFTUtils.hpp
- TEST_F
- TEST
- TEST_F
- Trajectory
- NeighborGraph
- ProgramOptions
- OutmolParser
- XdatcarHeader
- fuzz_utils.hpp
- BFSScratch
- Onetep File Parser
- Cell
- PDF Comparison Renderer
- Histogram Metadata
- PDF Histogram Renderer
- SvgHistogramRenderer
- TEST
- VASP XDATCAR Reader
- Main App Controller
- TEST_F
- TEST_F
- VaspParser
- Roboto
- renderComparisonSvg
- TEST
- TEST_F
- GromacsReader
- RDFCalculator.cpp
- TEST_F
- PYBIND11_MODULE
- renderComparisonPdf
- SvgComparisonRenderer
- TEST
- TEST_F
- TEST_F
- TEST
- TEST
- PlotController.cpp
- progress_callback
- Trajectory.cpp
- PlotConfig
- TEST
- TEST
- TEST
- TEST
- Constants.hpp
- TEST
- TEST
- GPUBond
- string
- StructureAnalyzer
- LocalEntropyCalculator.cpp
- TEST
- Cell.cpp
- QEReader
- GPUSQCalculator.cu
- QETrajectoryParser
- CNCalculator
- TEST_F
- TEST_F
- PhysicalData.hpp
- CP2KReader
- MappedFileFunctionalTests
- TEST_F
- TEST_F
- TEST
- TEST_F
- TEST
- VoronoiCalculator
- BaseReader
- CNACalculator.cpp
- TEST_F
- ArrowWriter.cpp
- BaseCalculator
- SteinhardtCalculator
- ReaderFactory
- DatasetWriteQuery
- PdfPlotter.hpp
- PyBaseCalculator
- AppController.cpp
- NiceScale
- TEST_F
- StructureFactorCalculator.cpp
- CellReader.cpp
- TEST_F
- GPUSQCalculator
- XRDCalculator::calculate
- .atomCount
- TEST_F
- WriterFactory
- GPUDistanceCalculator.cu
- HBondCalculator.cpp
- LammpsFrameParser
- DihedralCalculatorTests
- TEST
- DADCalculator
- DihedralCalculator
- LocalEntropyCalculator
- PADCalculator
- RDCalculator
- VDOSCalculator
- XRDCalculator
- BaseWriter
- TEST_F
- PdfPlotterTests
- SteinhardtCalculator.cpp
- GromacsReader.cpp
- TEST_F
- addFrame
- TEST_F
- RDFCalculator
- VACFCalculator
- renderTextAsPath
- CSVWriter
- Correlation Core Library
- UnionFind
- GPULattice
- PADCalculator.cpp
- ArcReader.cpp
- TEST_F
- wasm_bindings.cpp
- readTrajectory
- ClusterCalculator
- computeSingleAtomSteinhardt
- ArgBuilder
- TEST
- buildPartialsInfo
- XRDCalculator.cpp
- LammpsDumpReader::readTrajectory
- HyperuniformityCalculatorTests.cpp
- VoronoiCalculator::populateHistogram
- FileIOHandler
- CalculatorFactory
- HyperuniformityCalculator::calculate
- SDFCalculator
- MappedFile.hpp
- ArrowWriter
- HDF5Writer
- HistogramConfigs
- PartialInfo
- precomputePhases
- NeighborGraph.cpp
- MappedFileTests
- MockReader
- TEST_F
- AnalysisRunner
- HoverInfo
- AnalysisRunner.cpp
- GPUSearchGrid
- PrecomputedPhases
- ReciprocalBasis
- ThreadAccumulators
- TrajectoryTests.cpp
- File Reader Tests
- TEST
- FileWriter::write
- computeW6
- FileIOHandler.cpp
- GPUAtomData
- QVectorsData
- app.js
- Correlation: An Analysis Tool for Liquids and for Amorphous Solids
- hipLaunchKernelGGL
- main.cpp
- PlotController::requestPlotUpdate
- CSVWriter.cpp
- Coordinate Array Management
- CASTEP MD File Format
- BinningConfig
- SteinhardtParams
- PlotController::PlotController
- PlotSize
- ComparisonQuery
- __init__.py
- CliParserTests
- Common Neighbor Analysis (CNA)
- Steinhardt Parameters
- Code of Conduct
- Contributing Guidelines
- Roadmap to 4.0.0
- Energy State (E)
- Pressure (P)
- Temperature (T)
- CI Build Workflow
- Deploy Documentation Workflow
- Release Workflow
- Correlation Banner
- Demo Animation
- Radial Distribution Function Diagram
- Comparison Plots
- Correlation Logo
- RDF Plot Example
- correlation-analysis
- Correlation Build System
- Arrow Back Icon
- Arrow Drop Down Icon
- Arrow Drop Up Icon
- Arrow Left Icon
- Arrow Right Icon
- Calendar Today Icon
- Check Icon
- Chevron Backward Icon
- Chevron Forward Icon
- Close Icon
- Dark Mode Icon
- Edit Icon
- Keyboard Icon
- Light Mode Icon
- Menu Icon
- Pin Icon
- Remove Icon
- Save Icon
- Schedule Icon
- Constants.hpp
- Rule: Caveman Navigation Protocol (Token-Efficient Codebase Discovery)
- Expressions & Statements (ES.*)
- Error Handling (E.*)
- Templates & Generic Programming (T.*)
- Performance (Per.*)
- Philosophy & Interfaces (P.*, I.*)
- CliParser.hpp

## God Nodes (most connected - your core abstractions)
1. `Cell` - 221 edges
2. `DistributionFunctions` - 153 edges
3. `Trajectory` - 125 edges
4. `Histogram` - 71 edges
5. `TEST_F()` - 71 edges
6. `AppBackend` - 65 edges
7. `AnalysisSettings` - 53 edges
8. `StructureAnalyzer` - 53 edges
9. `NeighborGraph` - 51 edges
10. `PlotController` - 48 edges

## Surprising Connections (you probably didn't know these)
- `XRDCalculator::calculate()` --calls--> `calculatePartialIntegrands`  [INFERRED]
  src/calculators/XRDCalculator.cpp → include/calculators/XRDCalculator.hpp
- `ArcReader::read()` --calls--> `parseLine`  [INFERRED]
  src/readers/ArcReader.cpp → include/readers/ArcReader.hpp
- `readTrajectory()` --calls--> `isTrajectory`  [INFERRED]
  src/readers/FileReader.cpp → include/readers/BaseReader.hpp
- `GromacsReader::readTrajectory()` --calls--> `parseGroFrame`  [INFERRED]
  src/readers/GromacsReader.cpp → include/readers/GromacsReader.hpp
- `JOSS Logo` --references--> `Correlation: An Analysis Tool for Liquids and for Amorphous Solids`  [INFERRED]
  Images/joss-logo.png → paper/paper.md

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **CASTEP MD Test Data Suite** — tests_data_castep_md_clean, tests_data_castep_md_test, tests_fuzz_corpus_castep_md_minimal [INFERRED 0.90]
- **Material UI Icons** — ui_material_ui_icons_arrow_back, ui_material_ui_icons_arrow_drop_down, ui_material_ui_icons_arrow_drop_up, ui_material_ui_icons_arrow_left, ui_material_ui_icons_arrow_right, ui_material_ui_icons_calendar_today, ui_material_ui_icons_check, ui_material_ui_icons_chevron_backward, ui_material_ui_icons_chevron_forward, ui_material_ui_icons_close, ui_material_ui_icons_dark_mode, ui_material_ui_icons_edit, ui_material_ui_icons_keyboard, ui_material_ui_icons_light_mode, ui_material_ui_icons_menu, ui_material_ui_icons_pin, ui_material_ui_icons_remove, ui_material_ui_icons_save, ui_material_ui_icons_schedule [EXTRACTED 1.00]
- **CI/CD and Release Pipeline** — github_workflows_ci, github_workflows_docs, github_workflows_release [EXTRACTED 1.00]
- **Core Scientific Calculators** — concept_rdf, concept_steinhardt, concept_cna [EXTRACTED 1.00]
- **Cellulose Molecular Composition** — examples_cellulose_cellulose_md_carbon_atoms, examples_cellulose_cellulose_md_oxygen_atoms, examples_cellulose_cellulose_md_hydrogen_atoms [EXTRACTED 1.00]
- **Cellulose Molecular Structure** — examples_cellulose_cellulose_md_o_atoms, examples_cellulose_cellulose_md_h_atoms, examples_cellulose_cellulose_md_c_atoms [EXTRACTED 0.90]
- **Molecular Dynamics State** — examples_cellulose_cellulose_md_forces [INFERRED 0.85]
- **Cellulose Molecular Dynamics State** — examples_cellulose_cellulose_energy_state [INFERRED 0.90]
- **Cellulose Molecular Dynamics State** — examples_cellulose_cellulose_v_matrix, examples_cellulose_cellulose_energy [INFERRED 0.90]
- **Cellulose Molecular Dynamics Data** — examples_cellulose_cellulose_coordinates, examples_cellulose_cellulose_forces, examples_cellulose_cellulose_vibrational_modes [EXTRACTED 0.95]
- **Cellulose Computational Chemistry Data** — examples_cellulose_cellulose_coordinates, examples_cellulose_cellulose_vibrational_modes, examples_cellulose_cellulose_force_constants [EXTRACTED 0.95]
- **Cellulose Computational Chemistry Results** — examples_cellulose_cellulose_structure, examples_cellulose_cellulose_vibrational_modes, examples_cellulose_cellulose_force_constants, examples_cellulose_cellulose_energy_properties [EXTRACTED 0.95]
- **Core Analysis Architecture** — src_readers_obj, src_calculators_obj, src_writers_obj, src_correlation_lib [EXTRACTED 1.00]
- **Structural Distribution Functions** — concept_pdf, concept_pad, concept_rdf [EXTRACTED 1.00]

## Communities (256 total, 43 thin omitted)

### Community 0 - "PDF Generation Utilities"
Cohesion: 0.05
Nodes (118): pdf_object, FILE, determine_image_format(), dgets(), dstr_append(), dstr_append_data(), dstr_data(), dstr_ensure() (+110 more)

### Community 1 - "CLI Option Validation Tests"
Cohesion: 0.03
Nodes (67): AngleBinCannotExceed180, AngleBinMustBePositive, AngleBinOption, AtLeastOneOutputFormatMustBeEnabled, CrystallineMaterialDefaults, CsvEnableExplicit, CsvToggle, DefaultOutputBaseIsStemOfInput (+59 more)

### Community 2 - "ChiralityCalculator.cpp"
Cohesion: 0.14
Nodes (21): ChiralityCalculator, calculate, calculateFrame, computeSingleAtomChirality, string, accumulateHistogramMap(), addValueToHistogram(), BinningConfig (+13 more)

### Community 3 - "PresetManager.cpp"
Cohesion: 0.06
Nodes (51): AppController, AppWindow, vector, PresetController, handleDeletePreset, handleLoadPreset, handleMaterialTypeChanged, handleSavePreset (+43 more)

### Community 4 - "TrajectoryAnalyzer"
Cohesion: 0.06
Nodes (42): EndFrame, value, MaxFrames, value, StartFrame, value, real_t, vector (+34 more)

### Community 5 - "Atom"
Cohesion: 0.06
Nodes (40): AccessorsModifyStateCorrectly, AngleFunctionCalculatesNinetyDegrees, AngleFunctionClampsFloatingPointInaccuracies, AngleFunctionHandlesCoincidentAtoms, AngleFunctionHandlesCollinearAtoms, AngleFunctionHandlesNaNCoordinates, DistanceBetweenIdenticalAtomsIsZero, DistanceFunctionCalculatesCorrectly (+32 more)

### Community 6 - "KernelGenerationParams"
Cohesion: 0.20
Nodes (23): CompatibilityOverloadDerivesBinWidth, GenerateKernelNormalizesAndCalculatesCorrectly, fillBiweight(), fillBump(), fillCosine(), fillEpanechnikov(), fillGaussian(), fillTriweight() (+15 more)

### Community 7 - "XYZReader"
Cohesion: 0.07
Nodes (33): CommentData, string, vector, XYZReader, parseCommentLine, parseEnergy, parseLattice, parseProperties (+25 more)

### Community 8 - "LinearAlgebra.hpp"
Cohesion: 0.12
Nodes (22): CORRELATION_ALIGN, cross(), determinant(), distance(), dot(), size_t, T, U (+14 more)

### Community 9 - "vector"
Cohesion: 0.10
Nodes (9): atomic, map, string, vector, TrajectoryAnalyzer, writeHistogramToCSV, string, CSVWriter::writeAllCSVs() (+1 more)

### Community 10 - "AppBackend"
Cohesion: 0.06
Nodes (43): AppBackend, analysis_thread_func, calculateDynamicProperties, cancel_flag_, df_, getAshcroftWeights, getAtomCounts, getAvailableHistogramNames (+35 more)

### Community 11 - "SIMDUtils.hpp"
Cohesion: 0.07
Nodes (41): AtomRange, FactorialCorrectness, complex_exp_sum(), compute_dsq_block(), debye_sum(), dist_sq_scalar(), dot_block(), fill_position_block() (+33 more)

### Community 12 - "Cellulose Example Data"
Cohesion: 0.08
Nodes (41): Cellulose Example Data, Carbon Atoms (C), Hydrogen Atoms (H), Oxygen Atoms (O), Carbon 1, Cellulose Molecular Structure, Atomic Coordinates (R), Energy (E) (+33 more)

### Community 13 - "DistributionFunctions"
Cohesion: 0.06
Nodes (43): CalculateVACF_and_VDOS, CalculateVACF_GasLike, CalculateVACF_WithFrameRange, ComputeDiffusionCoefficientVACF_and_RelaxationTime, DistributionFunctionsDynamicProperties, DistributionFunctionsNonPhysicalOptions, DynamicsAnalyzerNonPhysicalInputs, DistributionFunctions (+35 more)

### Community 14 - "DistributionFunctions.cpp"
Cohesion: 0.10
Nodes (33): calculateAshcroftWeights, ensureNeighborsComputed, neighbors, smooth, KernelType, real_t, string, vector (+25 more)

### Community 15 - "ThreadLocalDistances"
Cohesion: 0.09
Nodes (38): DistanceTensor, real_t, vector, DistanceCalculator::compute(), FlatCellList, indices, offsets, SearchGridConfig (+30 more)

### Community 16 - "CliOptions"
Cohesion: 0.06
Nodes (32): CliOptions, angle_bin_width, csv, dihedral_bin_width, disable_groups, has_dihedral_bin, hdf5, hyper_samples (+24 more)

### Community 17 - "PlotController"
Cohesion: 0.05
Nodes (38): AppWindow, atomic, shared_ptr, size_t, string, thread, vector, PlotController (+30 more)

### Community 18 - "InputValidator.cpp"
Cohesion: 0.11
Nodes (34): AppErrors, AppController, AppWindow, InputValidator, backend_, controller_, handleCopyCliCommand, updateCliCommand (+26 more)

### Community 19 - "TEST_F"
Cohesion: 0.06
Nodes (35): CrossProductProperties, DeterminantAndInversion, DistanceFunction, InvertRoundTripGeneral, Matrix3AdditionSubtraction, Matrix3ArrayConversion, Matrix3ConstructorsAndAccessors, Matrix3Equality (+27 more)

### Community 20 - "AppBackend.cpp"
Cohesion: 0.13
Nodes (15): write, AppBackend::getAshcroftWeights(), AppBackend::getAvailableHistogramNames(), AppBackend::getBondCutoff(), AppBackend::getHistogram(), AppBackend::getRecommendedBondCutoffs(), AppBackend::getTimeStep(), AppBackend::load_file() (+7 more)

### Community 21 - "TEST_F"
Cohesion: 0.09
Nodes (22): AccessorsWork, AddAndScale, CalculateCoordinationNumber, CalculateRDF, ComputeMean, DefaultConstructorWorks, HandlesMissingPartialInAdd, add (+14 more)

### Community 22 - "CastepMdReader"
Cohesion: 0.11
Nodes (24): CastepMdReaderTests, CastepMdReader, parseAtomLine, parseEnergyLine, parseLatticeLine, read, readStructure, updateProgress (+16 more)

### Community 23 - "AnalysisSettings"
Cohesion: 0.06
Nodes (29): AnalysisSettings, active_calculators, angle_bin_width, cancel_flag, dihedral_bin_width, hyperuniformity_samples, lef_cutoff, lef_sigma (+21 more)

### Community 24 - "AppDefaults"
Cohesion: 0.06
Nodes (32): AppDefaults, ANGLE_BIN_WIDTH, ANGLE_BIN_WIDTH_CRYSTAL, ANGLE_BIN_WIDTH_LIQUID, LEF_CUTOFF, LEF_SIGMA, MSG_ANALYSIS_ABORTED, MSG_ANALYSIS_ENDED (+24 more)

### Community 25 - "CifReader.cpp"
Cohesion: 0.17
Nodes (30): ParseState, AsymmetricAtom, frac_pos, symbol, CifReader::read(), CifReader::readStructure(), CifReader::readTrajectory(), cleanCifValue() (+22 more)

### Community 26 - "TEST_F"
Cohesion: 0.06
Nodes (34): AddFrameAddsFrameToTrajectory, AddFrameThrowsOnAtomCountMismatch, AddFrameThrowsOnAtomOrderMismatch, AddFrameThrowsOnElementCountMismatch, AddFrameThrowsOnElementMismatch, CalculateVelocitiesComputesCorrectVelocities, CalculateVelocitiesDoesNotCrashOnEmptyTrajectory, CalculateVelocitiesHandlesPBC (+26 more)

### Community 27 - "FFTUtils.hpp"
Cohesion: 0.09
Nodes (27): AutocorrelateEmptyReturnsEmpty, AutocorrelateMatchesMathematicalDefinition, AutocorrelateReusesWorkspaceCorrectly, ComputeFFTHandlesEmptyInput, ComputeFFTHandlesNonPowerOfTwo, ComputeFFTHandlesPowerOfTwoAndInvert, ComputeFFTSizeOne, ComputeFFTThrowsOnNonPowerOfTwo (+19 more)

### Community 28 - "TEST_F"
Cohesion: 0.08
Nodes (26): CalculatePAD, EmptyCellThrows, EquilateralTriangle60, FullNormalizationCheck, Icosahedron_13Atoms, IcosahedronAnglesPAD, LinearGeometry180, MissingAnglesWhenCutoffIsTooSmall (+18 more)

### Community 29 - "TEST"
Cohesion: 0.08
Nodes (22): AcceptorOnly_NoHydrogens, BulkMetalNoHBonds, CoincidentDonorAndAcceptor_DoesNotCrash, CoincidentHAndDonor_DoesNotCrash, FluorineDonorAndAcceptor, HBondCalculatorTests, HBondCalculator, calculate (+14 more)

### Community 30 - "TEST_F"
Cohesion: 0.07
Nodes (29): AcosNumericalNoiseClamping, AddAtomRegistersNewElements, ConstructorThrowsOnZeroOrSingularVolume, ExtremelyLargeCell, ExtremelySmallCell, FindElementWorksCorrectly, HighAtomCount, inverse_lattice_vectors_ (+21 more)

### Community 31 - "Trajectory"
Cohesion: 0.05
Nodes (37): BasicUsage, CreateAnalyzerOutOfBoundsReturnsNullptr, computeMean, mutex, TrajectoryAnalyzer, FrameParser, mutex, real_t (+29 more)

### Community 32 - "NeighborGraph"
Cohesion: 0.10
Nodes (24): AddDirectedEdgeAndGetNeighbors, DenseAdjacencyMatrixMapping, DuplicateEdgesAreBothStored, vector, NeighborGraph, addDirectedEdge, adj_list_, getDenseAdjacencyMatrix (+16 more)

### Community 33 - "ProgramOptions"
Cohesion: 0.07
Nodes (29): KernelType, map, vector, ProgramOptions, active_calculators, angle_bin_width, bond_cutoffs_sq, dihedral_bin_width (+21 more)

### Community 34 - "OutmolParser"
Cohesion: 0.10
Nodes (21): function, ifstream, real_t, streampos, string, stringstream, vector, OutmolParser (+13 more)

### Community 35 - "XdatcarHeader"
Cohesion: 0.15
Nodes (22): function, real_t, shared_ptr, string, vector, extractLine(), findLineEnd(), parseXdatcarFrame() (+14 more)

### Community 36 - "fuzz_utils.hpp"
Cohesion: 0.11
Nodes (18): readTrajectory, readStructure, readTrajectory, readStructure, readTrajectory, readStructure, PresetManager::presetsDirectory(), LLVMFuzzerTestOneInput() (+10 more)

### Community 37 - "BFSScratch"
Cohesion: 0.14
Nodes (27): KingBFSSettings, PathEndpoints, RootSearchSettings, BFSScratch, cross_edges, dist, dist_king, local_cycles (+19 more)

### Community 38 - "Onetep File Parser"
Cohesion: 0.16
Nodes (17): function, ifstream, real_t, string, stringstream, OnetepDatParser, current_block_type, frac_flag (+9 more)

### Community 39 - "Cell"
Cohesion: 0.11
Nodes (20): BasicCalculation, IcosahedronAnglesDAD, Cell, energy_, getOrRegisterElement, lattice_parameters_, real_t, vector (+12 more)

### Community 40 - "PDF Comparison Renderer"
Cohesion: 0.09
Nodes (23): PdfComparisonRenderer, axis_col, bg_col, canvas_height, canvas_width, config, datasets, grid_col (+15 more)

### Community 41 - "Histogram Metadata"
Cohesion: 0.13
Nodes (16): Histogram, bins, compute_count, description, file_suffix, partials, smoothed_partials, title (+8 more)

### Community 42 - "PDF Histogram Renderer"
Cohesion: 0.08
Nodes (24): map, real_t, PdfHistogramRenderer, axis_col, bg_col, canvas_height, canvas_width, config (+16 more)

### Community 43 - "SvgHistogramRenderer"
Cohesion: 0.08
Nodes (27): HoverInfo, active, mouse_x, mouse_y, widget_height, widget_width, map, real_t (+19 more)

### Community 44 - "TEST"
Cohesion: 0.15
Nodes (12): BCC_Supercell_ProducesOutput, CNACalculatorTests, FCC_Supercell_ProducesNonEmptyResult, SingleIsolatedAtomReturnsEmptyHistogram, DeterministicResults, HistogramDimensionsAreConsistent, HistogramMetadataIsPopulated, NullNeighborsReturnsEmptyHistogram (+4 more)

### Community 45 - "VASP XDATCAR Reader"
Cohesion: 0.11
Nodes (19): FrameAtomCountConsistent, string, vector, XdatcarReader, readStructure, LatticeConsistentAcrossFrames, ParseThreeFrameTrajectory, PositionsDifferBetweenFrames (+11 more)

### Community 46 - "Main App Controller"
Cohesion: 0.08
Nodes (32): AnalysisRunner, AppController, analysis_runner_, file_io_handler_, getBondCutoffs, handleOptionsfromUI, handleOptionstoUI, input_validator_ (+24 more)

### Community 47 - "TEST_F"
Cohesion: 0.11
Nodes (17): string, vector, PdbReader, readStructure, readTrajectory, ReadFallbackElements, ReadSingleStructure, ReadTrajectoryFrames (+9 more)

### Community 48 - "TEST_F"
Cohesion: 0.11
Nodes (19): string, vector, VaspReader, read, readTrajectory, ParseCartesianCoordinates, ParseMultiSpecies, ParseSelectiveDynamics (+11 more)

### Community 49 - "VaspParser"
Cohesion: 0.18
Nodes (11): function, ifstream, pair, real_t, string, vector, VaspParser, file (+3 more)

### Community 50 - "Roboto"
Cohesion: 0.11
Nodes (22): Glyph, left, right, strokes, map, pair, vector, Roboto (+14 more)

### Community 51 - "renderComparisonSvg"
Cohesion: 0.27
Nodes (3): fmtScientific(), mapValue(), NearestPoint

### Community 52 - "TEST"
Cohesion: 0.11
Nodes (22): CalculatesMSDCorrectly, CalculatesVACFFromExampletraj, CalculatesVDOSCorrectly, ComputesDiffusionCoefficientMSD, ComputesDiffusionCoefficientVACF, ComputesRelaxationTime, DynamicsAnalyzerTests, HandlesEmptyAndInvalidTrajectories (+14 more)

### Community 53 - "TEST_F"
Cohesion: 0.13
Nodes (20): DefaultExecutesAllCalculators, DisableRadialAndScatteringGroups, DisableStructuralGroup, HelpReturnsZero, NoArgsReturnsNonZero, NonexistentFileReturnsNonZero, ShortHelpFlag, ShortVersionFlag (+12 more)

### Community 54 - "GromacsReader"
Cohesion: 0.12
Nodes (15): GromacsReader, parseGroFrame, readStructure, readTrajectory, string, vector, ReadMultiFrameTrajectory, ReadStructureReturnsLastFrame (+7 more)

### Community 55 - "RDFCalculator.cpp"
Cohesion: 0.19
Nodes (21): DistributionFunctions::calculateAshcroftWeights(), accumulateRawCounts(), map, real_t, string, getInversePartialKey(), getPartialKey(), normalizeDistributions() (+13 more)

### Community 56 - "TEST_F"
Cohesion: 0.08
Nodes (29): CalculatorInterfaceIsCorrect, HistogramMetadataIsCorrect, string, HyperuniformityCalculator, calculate, calculateFrame, LatticeHasLowerSlopeThanRandom, LatticeVarianceScalesAsR2 (+21 more)

### Community 57 - "PYBIND11_MODULE"
Cohesion: 0.12
Nodes (10): _correlation, mod, module_, init_core(), module_, init_io(), module_, init_math() (+2 more)

### Community 58 - "renderComparisonPdf"
Cohesion: 0.33
Nodes (6): parseHexColor(), renderComparisonPdf(), renderHistogramAsPdf(), pdf_add_circle(), pdf_add_line(), pdf_append_page()

### Community 59 - "SvgComparisonRenderer"
Cohesion: 0.08
Nodes (24): renderComparisonSvg(), renderHistogramAsSvg(), SvgComparisonRenderer, config, datasets, hover, kHeight, kWidth (+16 more)

### Community 60 - "TEST"
Cohesion: 0.12
Nodes (13): AtomsOutsideCutoff, ComputesPairwiseDistancesAndNeighborGraph, DistanceAcrossPeriodicBoundary, DistanceCalculatorTests, DistanceCalculator, calculateFrame, compute, string (+5 more)

### Community 61 - "TEST_F"
Cohesion: 0.20
Nodes (10): HexagonalClosePacked, lattice_vectors_, volume_, BodyCenteredCubic, FaceCenteredCubic, real_t, SimpleCubic, testing::Test (+2 more)

### Community 62 - "TEST_F"
Cohesion: 0.12
Nodes (16): CalculatesAndWritesSiliconDistributions, calculatePAD, calculateVACF, calculateVDOS, smoothAll, real_t, string, testing::Test (+8 more)

### Community 63 - "TEST"
Cohesion: 0.13
Nodes (15): CarReaderTests, CarReader, read, readStructure, readTrajectory, string, vector, ReadsStructure (+7 more)

### Community 64 - "TEST"
Cohesion: 0.13
Nodes (15): CellReaderTests, CellReader, read, readStructure, readTrajectory, string, vector, ReadsStructureLatticeAbc (+7 more)

### Community 65 - "PlotController.cpp"
Cohesion: 0.18
Nodes (19): executeSavePlot, requestPlotUpdate, RenderTaskData, SharedString, string, T, getComparisonKey(), PlotController::buildPlotConfigFromUI() (+11 more)

### Community 66 - "progress_callback"
Cohesion: 0.18
Nodes (19): setLatticeParameters, progress_callback, function, optional, real_t, string, parsePdbAtomLine(), parsePdbCrystLine() (+11 more)

### Community 67 - "Trajectory.cpp"
Cohesion: 0.14
Nodes (18): ensureMaterialized, getBondCutoffSQ, parser_, precomputeBondCutoffs, removeDuplicatedFrames, validateFrame, FrameParser, real_t (+10 more)

### Community 68 - "PlotConfig"
Cohesion: 0.10
Nodes (16): PlotConfig, fill_area, font_scale, height, line_width, palette, preset_size, show_grid (+8 more)

### Community 69 - "TEST"
Cohesion: 0.09
Nodes (29): ColumnLayout, string, vector, LammpsDumpReader, parseDumpFrame, readStructure, readTrajectory, LammpsDumpReaderTests (+21 more)

### Community 70 - "TEST"
Cohesion: 0.15
Nodes (13): string, vector, OnetepDatReader, read, readTrajectory, OnetepDatReaderTests, ReadsStructureCartesianAndBohr, ReadsStructureLatticeAbcAndPositionsFrac (+5 more)

### Community 71 - "TEST"
Cohesion: 0.13
Nodes (15): string, vector, OutmolReader, read, readStructure, readTrajectory, OutmolReaderTests, ReadsTrajectoryFormat1 (+7 more)

### Community 72 - "TEST"
Cohesion: 0.12
Nodes (16): ArcReaderTests, ArcReader, parseLine, read, readStructure, readTrajectory, updateProgress, string (+8 more)

### Community 73 - "Constants.hpp"
Cohesion: 0.12
Nodes (16): CalculatesCorrectAnglesForWater, CalculatesCorrectAngleWithPBC, CalculatesCorrectDihedralAngles, DistancesTensorIsCorrect, EnforcesNeighborSymmetry, FindsCorrectNeighborsForSilicon, FindsNeighborsBasedOnBondCutoff, FindsNoNeighborsForIsolatedAtom (+8 more)

### Community 74 - "TEST"
Cohesion: 0.14
Nodes (14): CifReaderTests, CifReader, read, readStructure, readTrajectory, string, vector, ReadsStructureWithSymmetryAndUncertainty (+6 more)

### Community 75 - "TEST"
Cohesion: 0.13
Nodes (12): ComputeDiffusionCoefficientMSD, ComputesCorrectMSDAndDeff, DynamicsAnalyzerMSDNonPhysicalInputs, FrameRangeSubset, string, MSDCalculator, calculate, calculateTrajectory (+4 more)

### Community 76 - "GPUBond"
Cohesion: 0.12
Nodes (19): __device__, __global__, distance_kernel(), GPUBinData, indices, offsets, GPUBond, distance (+11 more)

### Community 77 - "string"
Cohesion: 0.17
Nodes (10): pair, string, TextAnchor, tuple, vector, LabeledHistogram, hist, label (+2 more)

### Community 78 - "StructureAnalyzer"
Cohesion: 0.09
Nodes (20): AngleTensor, DihedralTensor, DistanceTensor, real_t, vector, StructureAnalyzer, angle_tensor_, bond_cutoffs_sq_ (+12 more)

### Community 79 - "LocalEntropyCalculator.cpp"
Cohesion: 0.27
Nodes (18): accumulateHistogramMap(), addValueToHistogram(), BinningConfig, d_val, max_val, min_val, computeSingleAtomEntropy(), copyPartialsToHistogram() (+10 more)

### Community 80 - "TEST"
Cohesion: 0.14
Nodes (11): AngleCalculatorTests, ComputesCorrect180DegreeAngle, ComputesCorrect60DegreeAngle, ComputesCorrect90DegreeAngle, AngleCalculator, calculateFrame, compute, string (+3 more)

### Community 81 - "Cell.cpp"
Cohesion: 0.18
Nodes (12): findElement, updateLattice, updateLatticeParametersFromVectors, Cell::Cell(), Cell::findElement(), Cell::getOrRegisterElement(), Cell::minimumImage(), Cell::setLatticeParameters() (+4 more)

### Community 82 - "QEReader"
Cohesion: 0.14
Nodes (12): string, vector, QEReader, readStructure, ReadsNonOrthogonalLattice, ReadsSingleFrame, string, testing::Test (+4 more)

### Community 83 - "GPUSQCalculator.cu"
Cohesion: 0.16
Nodes (16): averageBinnedSQ(), __global__, DeviceAtoms, x, y, z, DeviceQVectors, qx (+8 more)

### Community 84 - "QETrajectoryParser"
Cohesion: 0.16
Nodes (12): function, ifstream, string, vector, QEReader::readStructure(), QEReader::readTrajectory(), QETrajectoryParser, current_cell (+4 more)

### Community 85 - "CNCalculator"
Cohesion: 0.15
Nodes (9): CalculatesCorrectCoordinationNumbers, CNCalculatorTests, HighCoordinationFCC, CNCalculator, calculate, calculateFrame, string, IsolatedAtomHasZeroCoordination (+1 more)

### Community 86 - "TEST_F"
Cohesion: 0.07
Nodes (36): ComputeDsqBlockFloatMatchesScalar, ComputeDsqBlockMatchesScalar, DebyeSumDoubleAndFloat, DotBlockMatchesScalar, KahanSummationPrecision, mt19937, NormalizeRDFBinsFloatMatchesScalar, NormalizeRDFBinsMatchesScalar (+28 more)

### Community 87 - "TEST_F"
Cohesion: 0.14
Nodes (12): ConstructorInitializesCorrectly, HandlesBondCutoffsCorrectly, HandlesCalculatorToggleSignal, string, PopulatesRecommendedBondCutoffs, PopulatesTableAndDynamicProperties, SynchronizesOptionsToAndFromUI, AppControllerTests (+4 more)

### Community 88 - "PhysicalData.hpp"
Cohesion: 0.18
Nodes (14): GetAtomicFormFactorsCorrectly, GetAtomicMassCorrectly, GetCovalentRadiusCorrectly, ElementData, form_factors, mass, radius, symbol (+6 more)

### Community 89 - "CP2KReader"
Cohesion: 0.14
Nodes (12): CP2KReader, readStructure, readTrajectory, string, vector, CP2KReaderTests, data_dir_, ReadsSingleFrame (+4 more)

### Community 90 - "MappedFileFunctionalTests"
Cohesion: 0.12
Nodes (14): string, testing::Test, MappedFileFunctionalTests, content_a_, content_b_, file_a_path_, file_b_path_, test_dir_ (+6 more)

### Community 91 - "TEST_F"
Cohesion: 0.14
Nodes (14): BCC, EmptySystemOrNoNeighborsFillsPartialsWithZeros, FCC, HandlesAcosNumericalNoiseSafely, Icosahedral, SphericalHarmonics, map, real_t (+6 more)

### Community 92 - "TEST_F"
Cohesion: 0.14
Nodes (13): CalculateXRD, CalculateXRD_IntensityIsZeroAtThetaZero, CalculateXRD_InvalidBinWidth, CalculateXRD_InvalidInputsThrow, CalculateXRD_ThrowsIfNoRDF, CalculateXRDCubicCell, calculateRDF, calculateXRD (+5 more)

### Community 93 - "TEST"
Cohesion: 0.17
Nodes (9): CalculatorFactoryTests, GetRegisteredCalculators, LookupStandardCalculators, RegisterAndLookupCustomCalculator, LookupNonExistentReturnsNullptr, SingletonInstanceIsUnique, string, MockCalculator (+1 more)

### Community 94 - "TEST_F"
Cohesion: 0.10
Nodes (20): DetermineFileTypeExtensionlessVasp, ReadArcFileCorrectly, ReadArcFileDuplicatedFrames, ReadCarFileCorrectly, ReadCastepMdCorrectly, ReadCellFileCorrectly, ReadCelluloseExample, ReadCifFileCorrectly (+12 more)

### Community 95 - "TEST"
Cohesion: 0.15
Nodes (12): GetRegisteredWriters, GetWriterByName, GetWriterForExtension, RegisterAndLookupCustomWriter, LookupEmptyExtensionReturnsNullptr, LookupNonExistentReturnsNullptr, SingletonInstanceIsUnique, string (+4 more)

### Community 96 - "VoronoiCalculator"
Cohesion: 0.17
Nodes (10): string, VoronoiCalculator, buildSignatureMap, calculate, calculateFrame, computeVoronoiCells, makeHistogram, populateHistogram (+2 more)

### Community 97 - "BaseReader"
Cohesion: 0.16
Nodes (9): BaseReader, getExtensions, getName, isTrajectory, readStructure, readTrajectory, string, vector (+1 more)

### Community 98 - "CNACalculator.cpp"
Cohesion: 0.24
Nodes (15): set, buildCNAHistogram(), buildCommonNeighborAdjacency(), CNACalculator::calculate(), CNACalculator::calculateFrame(), CommonNeighborAdjacency, adjacency_list, bond_count (+7 more)

### Community 99 - "TEST_F"
Cohesion: 0.15
Nodes (13): DetectsSingleSquare, DetectsSingleTriangle, EmptyGraphReturnsNoRings, MotifFinder, extractCycles, findRings, IsolatedNodesReturnsNoRings, MaxRingSizeExcludesLargerRings (+5 more)

### Community 100 - "ArrowWriter.cpp"
Cohesion: 0.23
Nodes (14): FieldVector, writeHistogramToParquet, addFloatColumn(), ArrowWriter::writeAllParquet(), ArrowWriter::writeHistogramToParquet(), Array, map, real_t (+6 more)

### Community 101 - "BaseCalculator"
Cohesion: 0.14
Nodes (11): BaseCalculator, getDescription, getGroup, getName, getShortName, isFrameCalculator, isTrajectoryCalculator, CalculatorFactory::getCalculator() (+3 more)

### Community 102 - "SteinhardtCalculator"
Cohesion: 0.21
Nodes (5): string, SteinhardtCalculator, calculate, calculateFrame, wigner3j

### Community 103 - "ReaderFactory"
Cohesion: 0.20
Nodes (9): map, unique_ptr, vector, ReaderFactory, extension_map_, getAllExtensions, getReaderForExtension, readers_ (+1 more)

### Community 104 - "DatasetWriteQuery"
Cohesion: 0.21
Nodes (13): Group, File, string, vector, DatasetWriteQuery, bin_ds_name, bin_unit, data_unit (+5 more)

### Community 105 - "PdfPlotter.hpp"
Cohesion: 0.16
Nodes (12): blendColor(), BlendParams, bg, fg, drawPdfText(), fmtScientificPdf(), TextAnchor, PdfPoint (+4 more)

### Community 106 - "PyBaseCalculator"
Cohesion: 0.19
Nodes (4): module_, string, init_calculators(), PyBaseCalculator

### Community 107 - "AppController.cpp"
Cohesion: 0.11
Nodes (13): CalculateSF_EmptyCellThrows, CalculateSF_InvalidInputsThrow, CalculatesSimpleCubicBraggPeak, DimerProducesValidSQ, HomonuclearClusterPartialsPresent, string, StructureFactorCalculator, calculateFrame (+5 more)

### Community 108 - "NiceScale"
Cohesion: 0.18
Nodes (8): DataRange, max, min, NiceScale, max, min, spacing, ticks

### Community 109 - "TEST_F"
Cohesion: 0.17
Nodes (10): LoadAllSorting, MalformedJsonHandling, MissingKeysFallback, PresetNameSanitization, RoundTripSerialization, SaveAndLoadAll, SerializationWithSpecialCharacters, testing::Test (+2 more)

### Community 110 - "StructureFactorCalculator.cpp"
Cohesion: 0.05
Nodes (68): LocalAccumulatorTuple, buildPartialsInfo(), buildTypeBlocks(), combineAccumulators(), computeReciprocalBasis(), CoordinateArrays, x, y (+60 more)

### Community 111 - "CellReader.cpp"
Cohesion: 0.31
Nodes (13): array(), CellReader::read(), CellReader::readStructure(), CellReader::readTrajectory(), function, real_t, string, stringstream (+5 more)

### Community 112 - "TEST_F"
Cohesion: 0.22
Nodes (8): CelluloseRingDistribution, ComputeMotif, neighbor_graph_, InvalidMaxRingSize, testing::Test, RDTests, graph, TEST_F()

### Community 113 - "GPUSQCalculator"
Cohesion: 0.21
Nodes (4): GPUSQCalculator, calculateFrame, has_gpu_, string

### Community 114 - "XRDCalculator::calculate"
Cohesion: 0.21
Nodes (12): BinWidth, value, real_t, MaxTheta, value, MinTheta, value, Wavelength (+4 more)

### Community 115 - ".atomCount"
Cohesion: 0.26
Nodes (10): Cp2kParserState, has_box, parsing_cell, parsing_coords, CP2KReader::readStructure(), CP2KReader::readTrajectory(), function, string (+2 more)

### Community 116 - "TEST_F"
Cohesion: 0.10
Nodes (20): AccumulationPrecisionLoss, Matrix3Determinant, Matrix3Inversion, testing::Test, PrecisionTests, d_x1, d_x2, d_y1 (+12 more)

### Community 117 - "WriterFactory"
Cohesion: 0.17
Nodes (11): map, string, unique_ptr, vector, WriterFactory, extension_map_, getWriter, getWriterForExtension (+3 more)

### Community 118 - "GPUDistanceCalculator.cu"
Cohesion: 0.30
Nodes (10): compute_distances_gpu(), DistanceTensor, real_t, vector, flatten_bond_cutoffs(), GPUPosition, x, y (+2 more)

### Community 119 - "HBondCalculator.cpp"
Cohesion: 0.26
Nodes (12): checkAcceptorsForHydrogen(), real_t, string, vector, findBondedHydrogens(), findHydrogenBonds(), HBondCalculator::calculate(), HBondCalculator::calculateFrame() (+4 more)

### Community 120 - "LammpsFrameParser"
Cohesion: 0.15
Nodes (15): BuildBCCLatticeAndVerifyDensity, BuildFCCLatticeAndVerifyPBCDistances, atoms_, elements_, minimumImage, wrapPositions, getFrame, AppBackend::getAtomCounts() (+7 more)

### Community 121 - "DihedralCalculatorTests"
Cohesion: 0.20
Nodes (9): ComputesCorrect0DegreeDihedral, ComputesCorrect180DegreeDihedral, ComputesCorrect90DegreeDihedral, HandlesCoincidentCentralBondSafely, testing::Test, DihedralCalculatorTests, cell, graph (+1 more)

### Community 122 - "TEST"
Cohesion: 0.18
Nodes (11): GetAllExtensions, GetReaderForExtension, GetRegisteredReaders, ReaderFactoryTests, RegisterAndLookupCustomReader, SniffsCP2KFromOutFile, SniffsQuantumEspressoFromOutFile, LookupEmptyExtensionReturnsNullptr (+3 more)

### Community 123 - "DADCalculator"
Cohesion: 0.24
Nodes (4): DADCalculator, calculate, calculateFrame, string

### Community 124 - "DihedralCalculator"
Cohesion: 0.24
Nodes (4): DihedralCalculator, calculateFrame, compute, string

### Community 125 - "LocalEntropyCalculator"
Cohesion: 0.24
Nodes (4): string, LocalEntropyCalculator, calculate, calculateFrame

### Community 126 - "PADCalculator"
Cohesion: 0.24
Nodes (4): string, PADCalculator, calculate, calculateFrame

### Community 127 - "RDCalculator"
Cohesion: 0.24
Nodes (4): string, RDCalculator, calculate, calculateFrame

### Community 128 - "VDOSCalculator"
Cohesion: 0.16
Nodes (8): real_t, string, VDOSCalculator, calculate, calculateTrajectory, VDOSParams, max_imag_freq, max_real_freq

### Community 129 - "XRDCalculator"
Cohesion: 0.24
Nodes (5): string, XRDCalculator, calculate, calculateFrame, calculatePartialIntegrands

### Community 130 - "BaseWriter"
Cohesion: 0.24
Nodes (9): BaseWriter, getExtensions, getName, write, string, unique_ptr, WriterFactory::getWriter(), WriterFactory::getWriterForExtension() (+1 more)

### Community 131 - "TEST_F"
Cohesion: 0.21
Nodes (13): AngleCalculator::calculateFrame(), AngleCalculator::compute(), AngleScratch, dots, nb_dist, nb_x, nb_y, nb_z (+5 more)

### Community 132 - "PdfPlotterTests"
Cohesion: 0.24
Nodes (8): RendersComparisonPdfCorrectly, RendersEmptyHistogramAsPdfGracefully, RendersValidHistogramAsPdfCorrectly, string, testing::Test, PdfPlotterTests, temp_pdf_path, TEST_F()

### Community 133 - "SteinhardtCalculator.cpp"
Cohesion: 0.55
Nodes (10): accumulateHistogramMap(), addValueToHistogram(), copyPartialsToHistogram(), map, string, vector, initHistogramMap(), normalizeHistogramMap() (+2 more)

### Community 134 - "GromacsReader.cpp"
Cohesion: 0.40
Nodes (10): function, string, extractLine(), findLineEnd(), GromacsReader::parseGroFrame(), GromacsReader::readStructure(), GromacsReader::readTrajectory(), scanNextFrame() (+2 more)

### Community 135 - "TEST_F"
Cohesion: 0.22
Nodes (8): AchiralCoplanarMotif, ChiralLeftHandedMotif, ChiralRightHandedMotif, HistogramDistribution, ChiralityCalculatorTests, real_t, testing::Test, TEST_F()

### Community 136 - "addFrame"
Cohesion: 0.21
Nodes (14): calculateDihedralAngle(), DihedralTensor, enumerable_thread_specific, optional, real_t, vector, DihedralCalculator::calculateFrame(), DihedralCalculator::compute() (+6 more)

### Community 137 - "TEST_F"
Cohesion: 0.20
Nodes (8): DataUriLoadingFailsAsExpected, FailsWithDotExtension, HandlesInvalidSvgFromEmbeddedData, LoadsBasicSvgFromEmbeddedData, LoadsSvgWithGradientFromEmbeddedData, testing::Test, SlintSvgIntegrationTests, TEST_F()

### Community 138 - "RDFCalculator"
Cohesion: 0.27
Nodes (4): string, RDFCalculator, calculate, calculateFrame

### Community 139 - "VACFCalculator"
Cohesion: 0.27
Nodes (4): string, VACFCalculator, calculate, calculateTrajectory

### Community 140 - "renderTextAsPath"
Cohesion: 0.19
Nodes (11): string, ReaderExtensionQuery, extension, filename, string, unique_ptr, vector, ReaderFactory::getAllExtensions() (+3 more)

### Community 141 - "CSVWriter"
Cohesion: 0.23
Nodes (6): CSVWriter, writeAllCSVs, string, vector, module_, init_writers()

### Community 142 - "Correlation Core Library"
Cohesion: 0.20
Nodes (10): Python Bindings (_correlation), Calculators Library, Correlation CLI, Correlation GUI, Correlation Core Library, Readers Library, Writers Library, Fuzz Tests (+2 more)

### Community 143 - "UnionFind"
Cohesion: 0.31
Nodes (5): ClusterCalculator::calculateFrame(), vector, UnionFind, parent, sz

### Community 144 - "GPULattice"
Cohesion: 0.20
Nodes (10): GPULattice, v0_x, v0_y, v0_z, v1_x, v1_y, v1_z, v2_x (+2 more)

### Community 145 - "PADCalculator.cpp"
Cohesion: 0.29
Nodes (9): real_t, vector, PADCalculator::calculate(), PADCalculator::calculateFrame(), PADSettings, bin_width, num_bins, theta_cut (+1 more)

### Community 146 - "ArcReader.cpp"
Cohesion: 0.36
Nodes (9): ArcReader::parseLine(), ArcReader::read(), ArcReader::readStructure(), ArcReader::readTrajectory(), ArcReader::updateProgress(), function, streampos, string (+1 more)

### Community 147 - "TEST_F"
Cohesion: 0.17
Nodes (11): All Rules, C++ Coding Standards (C++ Core Guidelines), Constants & Immutability (Con.*), Cross-Cutting Principles, Enumerations (Enum.*), Key Rules, Key Rules, Quick Reference Checklist (+3 more)

### Community 148 - "wasm_bindings.cpp"
Cohesion: 0.33
Nodes (8): correlation_wasm, string, EMSCRIPTEN_BINDINGS(), getBinsJS(), getPartialJS(), getPartialKeysJS(), readFromBuffer(), val

### Community 149 - "readTrajectory"
Cohesion: 0.56
Nodes (8): FileType, function, string, determineFileType(), findReaderForFile(), findReaderForType(), readStructure(), readTrajectory()

### Community 151 - "computeSingleAtomSteinhardt"
Cohesion: 0.22
Nodes (9): sphericalHarmonic, computeSingleAtomSteinhardt(), GlobalSteinhardtFactors, global_Q4_factor, global_Q6_factor, SingleAtomSteinhardt, Q4, Q6 (+1 more)

### Community 152 - "ArgBuilder"
Cohesion: 0.28
Nodes (6): initializer_list, ArgBuilder, argv, storage, string, vector

### Community 153 - "TEST"
Cohesion: 0.22
Nodes (9): RendersComparisonOverlayCorrectly, RendersDarkThemeCorrectly, RendersEmptyHistogramGracefully, RendersShadedCurveCorrectly, RendersValidHistogramCorrectly, RendersWithHover2DNearestSnapping, RendersWithHoverActive, SvgPlotterTests (+1 more)

### Community 154 - "buildPartialsInfo"
Cohesion: 0.18
Nodes (10): 1. Minimal Code Modifications, 2. Concise Diagnostics, 3. Structural Scannability, Behavioral Blueprint, Caveman Communication & Token Economy Protocol, DO NOT USE (Filler Phrases), DO USE (Direct Answers), Limitations (+2 more)

### Community 155 - "XRDCalculator.cpp"
Cohesion: 0.36
Nodes (8): map, real_t, string, vector, XRDCalculator::calculateConcentrations(), XRDCalculator::calculateFrame(), XRDCalculator::calculatePartialIntegrands(), XRDCalculator::getAtomicFormFactor()

### Community 156 - "LammpsDumpReader::readTrajectory"
Cohesion: 0.24
Nodes (4): CNACalculator, calculate, calculateFrame, string

### Community 157 - "HyperuniformityCalculatorTests.cpp"
Cohesion: 0.31
Nodes (3): color(), size_t, Palette

### Community 158 - "VoronoiCalculator::populateHistogram"
Cohesion: 0.21
Nodes (12): BinRange, CellData, pair, real_t, string, vector, processCellTopology(), VoronoiCalculator::buildSignatureMap() (+4 more)

### Community 159 - "FileIOHandler"
Cohesion: 0.17
Nodes (14): AnalysisRunner, analysis_thread_, handleRunAnalysis, updateProgress, AppController, AppWindow, thread, AppController (+6 more)

### Community 160 - "CalculatorFactory"
Cohesion: 0.25
Nodes (7): CalculatorFactory, calculators_, getCalculator, registerCalculator, unique_ptr, vector, registerTypeSafe()

### Community 161 - "HyperuniformityCalculator::calculate"
Cohesion: 0.14
Nodes (11): real_t, HyperuniformityParams, num_samples, r_bin_width, map, mt19937_64, real_t, string (+3 more)

### Community 163 - "MappedFile.hpp"
Cohesion: 0.09
Nodes (17): EnforceSizeLimitCheck, size_t, MappedFile, release(), size(), MapsValidFileSuccessfully, MoveAssignmentOperatorTransfersOwnership, MoveConstructorTransfersOwnership (+9 more)

### Community 164 - "ArrowWriter"
Cohesion: 0.36
Nodes (4): ArrowWriter, writeAllParquet, string, vector

### Community 165 - "HDF5Writer"
Cohesion: 0.36
Nodes (4): HDF5Writer, writeHDF, string, vector

### Community 166 - "HistogramConfigs"
Cohesion: 0.25
Nodes (8): HistogramConfigs, bins_Q, bins_W, dQ, dW, Q_max, W_max, W_min

### Community 167 - "PartialInfo"
Cohesion: 0.25
Nodes (7): Approach, Do not use this skill when, Focus Areas, Instructions, Limitations, Output, Use this skill when

### Community 168 - "precomputePhases"
Cohesion: 0.46
Nodes (7): real_t, vector, DADCalculator::calculate(), DADCalculator::calculateFrame(), initializeHistogram(), normalizeAndScale(), processDihedralAngles()

### Community 169 - "NeighborGraph.cpp"
Cohesion: 0.12
Nodes (13): AtomIndex, id, AtomID, real_t, Neighbor, distance, index, r_ij (+5 more)

### Community 170 - "MappedFileTests"
Cohesion: 0.29
Nodes (6): 1. Static Analysis & Formatting (Clang Tooling), 2. Resource & Memory Safety (RAII & Smart Pointers), 3. Modern C++ Architecture (C++20/C++23), 4. Doxygen Documentation Standard, Example, Rule: C++ Code Style, Memory Safety, and Modern Standards

### Community 171 - "MockReader"
Cohesion: 0.29
Nodes (6): BasicClustering, ClusterCalculatorTests, EmptyCell, calculateFrame, SingleGiantCluster, TEST()

### Community 172 - "TEST_F"
Cohesion: 0.33
Nodes (6): Anti-Patterns, Class Hierarchy, Classes & Class Hierarchies (C.*), Key Rules, Rule of Five, Rule of Zero

### Community 173 - "AnalysisRunner"
Cohesion: 0.33
Nodes (5): Execution Pipeline, Phase 1: Sanitizer Flag Injection, Phase 2: Runtime Execution & Trap Capture, Phase 3: Root-Cause Remediation & Resolution Loop, Sanitizer Validator

### Community 174 - "HoverInfo"
Cohesion: 0.53
Nodes (5): App, applyMaterialDefaults(), parseArgs(), printUsage(), validateOptions()

### Community 175 - "AnalysisRunner.cpp"
Cohesion: 0.38
Nodes (5): AnalysisRunner::AnalysisRunner(), AnalysisRunner::updateProgress(), AppController, AppWindow, string

### Community 176 - "GPUSearchGrid"
Cohesion: 0.29
Nodes (7): GPUSearchGrid, K_x, K_y, K_z, max_dx, max_dy, max_dz

### Community 177 - "PrecomputedPhases"
Cohesion: 0.53
Nodes (5): CarReader::read(), CarReader::readStructure(), CarReader::readTrajectory(), function, string

### Community 178 - "ReciprocalBasis"
Cohesion: 0.40
Nodes (4): 1. Workspace Context, 2. Rule Hierarchy & Discovery Strategy, 3. Delegation & Skill Invocation, Global Agent Directives

### Community 179 - "ThreadAccumulators"
Cohesion: 0.40
Nodes (5): Anti-Patterns, Functions (F.*), Key Rules, Parameter Passing, Pure Functions and constexpr

### Community 180 - "TrajectoryTests.cpp"
Cohesion: 0.40
Nodes (5): Anti-Patterns, Key Rules, RAII Pattern, Resource Management (R.*), Smart Pointer Usage

### Community 181 - "File Reader Tests"
Cohesion: 0.40
Nodes (5): Anti-Patterns, Concurrency & Parallelism (CP.*), Key Rules, Multiple Mutexes, Safe Locking

### Community 182 - "TEST"
Cohesion: 0.33
Nodes (5): CalculateSDF, calculateFrame, MultiAtomSDFHasNonzeroDensity, SDFCalculatorTests, TEST()

### Community 183 - "FileWriter::write"
Cohesion: 0.40
Nodes (5): writeSummaryFile, string, FileWriter::FileWriter(), FileWriter::write(), FileWriter::writeSummaryFile()

### Community 184 - "computeW6"
Cohesion: 0.29
Nodes (8): SphericalAngles, computeW6(), complex, real_t, SteinhardtCalculator::sphericalHarmonic(), SteinhardtCalculator::wigner3j(), Wigner6Table, table

### Community 185 - "FileIOHandler.cpp"
Cohesion: 0.47
Nodes (3): AppController, AppWindow, FileIOHandler::FileIOHandler()

### Community 186 - "GPUAtomData"
Cohesion: 0.33
Nodes (6): GPUAtomData, atom_bin, element_ids, wrapped_x, wrapped_y, wrapped_z

### Community 187 - "QVectorsData"
Cohesion: 0.33
Nodes (6): vector, QVectorsData, qmag, qx, qy, qz

### Community 188 - "app.js"
Cohesion: 0.47
Nodes (4): COLORS, drawChart(), renderPlot(), runAnalysis()

### Community 189 - "Correlation: An Analysis Tool for Liquids and for Amorphous Solids"
Cohesion: 0.40
Nodes (5): Plane Angle Distribution (g(theta)), Pair Distribution Function (g(r)), Radial Distribution Function (RDF), JOSS Logo, Correlation: An Analysis Tool for Liquids and for Amorphous Solids

### Community 190 - "hipLaunchKernelGGL"
Cohesion: 0.40
Nodes (4): cudaStream_t, dim3, hipLaunchKernelGGL(), K

### Community 191 - "main.cpp"
Cohesion: 0.50
Nodes (4): HINSTANCE, LPSTR, main(), WinMain()

### Community 192 - "PlotController::requestPlotUpdate"
Cohesion: 0.40
Nodes (5): buildPlotConfigFromUI, executePlotRender, isPlotCacheHit, updateTableData, PlotController::requestPlotUpdate()

### Community 193 - "CSVWriter.cpp"
Cohesion: 0.40
Nodes (5): Anti-Patterns, Header Guard, Key Rules, Naming Conventions, Source Files & Naming (SF.*, NL.*)

### Community 194 - "Coordinate Array Management"
Cohesion: 0.40
Nodes (4): Phase 1: Environment & Build Detection, Phase 2: The Self-Correction Compilation Loop, Phase 3: Test Execution & Verification, Test-and-Tea Loop

### Community 195 - "CASTEP MD File Format"
Cohesion: 0.50
Nodes (4): CASTEP MD File Format, CASTEP MD Clean Data, CASTEP MD Test Data, CASTEP MD Minimal Fuzz Corpus

### Community 196 - "BinningConfig"
Cohesion: 0.50
Nodes (4): BinningConfig, d_val, max_val, min_val

### Community 197 - "SteinhardtParams"
Cohesion: 0.50
Nodes (4): SteinhardtParams, Q4, Q6, W6_hat

### Community 198 - "PlotController::PlotController"
Cohesion: 0.67
Nodes (3): handleUpdateTimer, AppWindow, PlotController::PlotController()

### Community 199 - "PlotSize"
Cohesion: 0.67
Nodes (3): PlotSize, height, width

### Community 200 - "ComparisonQuery"
Cohesion: 0.31
Nodes (6): ComparisonQuery, filepath, key, pair, string, vector

### Community 248 - "Constants.hpp"
Cohesion: 0.40
Nodes (3): real_t, VDOSCalculator::calculate(), VDOSCalculator::calculateTrajectory()

### Community 249 - "Rule: Caveman Navigation Protocol (Token-Efficient Codebase Discovery)"
Cohesion: 0.50
Nodes (3): 1. Core Directives, 2. Execution Workflow, Rule: Caveman Navigation Protocol (Token-Efficient Codebase Discovery)

### Community 250 - "Expressions & Statements (ES.*)"
Cohesion: 0.50
Nodes (4): Anti-Patterns, Expressions & Statements (ES.*), Initialization, Key Rules

### Community 251 - "Error Handling (E.*)"
Cohesion: 0.50
Nodes (4): Anti-Patterns, Error Handling (E.*), Exception Hierarchy, Key Rules

### Community 252 - "Templates & Generic Programming (T.*)"
Cohesion: 0.50
Nodes (4): Anti-Patterns, Concepts (C++20), Key Rules, Templates & Generic Programming (T.*)

### Community 253 - "Performance (Per.*)"
Cohesion: 0.50
Nodes (4): Anti-Patterns, Guidelines, Key Rules, Performance (Per.*)

### Community 254 - "Philosophy & Interfaces (P.*, I.*)"
Cohesion: 0.50
Nodes (4): DO, DON'T, Key Rules, Philosophy & Interfaces (P.*, I.*)

## Knowledge Gaps
- **833 isolated node(s):** `value`, `value`, `value`, `TrajectoryAnalyzer`, `r_max` (+828 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **43 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `Cell` connect `Cell` to `ChiralityCalculator.cpp`, `Atom`, `XYZReader`, `LinearAlgebra.hpp`, `vector`, `AppBackend`, `DistributionFunctions`, `DistributionFunctions.cpp`, `ThreadLocalDistances`, `TEST_F`, `CastepMdReader`, `CifReader.cpp`, `TEST_F`, `TEST_F`, `TEST`, `TEST_F`, `Trajectory`, `NeighborGraph`, `OutmolParser`, `XdatcarHeader`, `Onetep File Parser`, `Histogram Metadata`, `TEST`, `VASP XDATCAR Reader`, `TEST_F`, `TEST_F`, `VaspParser`, `TEST`, `GromacsReader`, `RDFCalculator.cpp`, `TEST_F`, `PYBIND11_MODULE`, `TEST`, `TEST_F`, `TEST_F`, `progress_callback`, `Trajectory.cpp`, `TEST`, `Constants.hpp`, `TEST`, `StructureAnalyzer`, `LocalEntropyCalculator.cpp`, `TEST`, `Cell.cpp`, `GPUSQCalculator.cu`, `QETrajectoryParser`, `CNCalculator`, `TEST_F`, `TEST_F`, `TEST_F`, `TEST_F`, `VoronoiCalculator`, `BaseReader`, `CNACalculator.cpp`, `TEST_F`, `AppController.cpp`, `StructureFactorCalculator.cpp`, `CellReader.cpp`, `XRDCalculator::calculate`, `.atomCount`, `GPUDistanceCalculator.cu`, `HBondCalculator.cpp`, `LammpsFrameParser`, `DihedralCalculatorTests`, `DihedralCalculator`, `TEST_F`, `SteinhardtCalculator.cpp`, `GromacsReader.cpp`, `TEST_F`, `addFrame`, `PADCalculator.cpp`, `ArcReader.cpp`, `wasm_bindings.cpp`, `readTrajectory`, `XRDCalculator.cpp`, `VoronoiCalculator::populateHistogram`, `HyperuniformityCalculator::calculate`, `precomputePhases`, `MockReader`, `PrecomputedPhases`, `TEST`?**
  _High betweenness centrality (0.256) - this node is a cross-community bridge._
- **Why does `DistributionFunctions` connect `DistributionFunctions` to `VDOSCalculator`, `ChiralityCalculator.cpp`, `TEST_F`, `addFrame`, `vector`, `AppBackend`, `CSVWriter`, `DistributionFunctions.cpp`, `UnionFind`, `ThreadLocalDistances`, `PADCalculator.cpp`, `wasm_bindings.cpp`, `TEST_F`, `AnalysisSettings`, `XRDCalculator.cpp`, `TEST_F`, `Trajectory`, `HyperuniformityCalculator::calculate`, `ArrowWriter`, `HDF5Writer`, `Cell`, `precomputePhases`, `Histogram Metadata`, `MockReader`, `TEST`, `FileWriter::write`, `TEST_F`, `TEST_F`, `string`, `StructureAnalyzer`, `GPUSQCalculator.cu`, `CNCalculator`, `TEST_F`, `TEST`, `CNACalculator.cpp`, `ArrowWriter.cpp`, `BaseCalculator`, `SteinhardtCalculator`, `DatasetWriteQuery`, `PdfPlotter.hpp`, `PyBaseCalculator`, `AppController.cpp`, `StructureFactorCalculator.cpp`, `GPUSQCalculator`, `HBondCalculator.cpp`, `Constants.hpp`, `DADCalculator`, `DihedralCalculator`, `LocalEntropyCalculator`, `PADCalculator`, `RDCalculator`?**
  _High betweenness centrality (0.188) - this node is a cross-community bridge._
- **Why does `Trajectory` connect `Trajectory` to `TrajectoryAnalyzer`, `GromacsReader.cpp`, `XYZReader`, `vector`, `AppBackend`, `DistributionFunctions`, `DistributionFunctions.cpp`, `ArcReader.cpp`, `wasm_bindings.cpp`, `readTrajectory`, `CastepMdReader`, `AnalysisSettings`, `TEST_F`, `CifReader.cpp`, `TEST_F`, `TEST_F`, `OutmolParser`, `XdatcarHeader`, `Onetep File Parser`, `Cell`, `VASP XDATCAR Reader`, `TEST_F`, `PrecomputedPhases`, `VaspParser`, `TEST`, `GromacsReader`, `PYBIND11_MODULE`, `TEST_F`, `progress_callback`, `Trajectory.cpp`, `TEST`, `Constants.hpp`, `TEST`, `QETrajectoryParser`, `PhysicalData.hpp`, `LammpsFrameParser`, `TEST_F`, `TEST_F`, `BaseReader`, `BaseCalculator`, `PyBaseCalculator`, `CellReader.cpp`, `TEST_F`, `.atomCount`, `Constants.hpp`?**
  _High betweenness centrality (0.096) - this node is a cross-community bridge._
- **Are the 29 inferred relationships involving `Cell` (e.g. with `CP2KReader::readTrajectory()` and `readTrajectory()`) actually correct?**
  _`Cell` has 29 INFERRED edges - model-reasoned connections that need verification._
- **Are the 4 inferred relationships involving `DistributionFunctions` (e.g. with `TEST_F()` and `TEST_F()`) actually correct?**
  _`DistributionFunctions` has 4 INFERRED edges - model-reasoned connections that need verification._
- **Are the 7 inferred relationships involving `Trajectory` (e.g. with `TEST()` and `TEST_F()`) actually correct?**
  _`Trajectory` has 7 INFERRED edges - model-reasoned connections that need verification._
- **What connects `value`, `value`, `value` to the rest of the system?**
  _833 weakly-connected nodes found - possible documentation gaps or missing edges._