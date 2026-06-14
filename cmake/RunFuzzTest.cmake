# cmake/RunFuzzTest.cmake
# Wraps execution of a fuzz target to read from the seed corpus directory,
# write temporary coverage files to a temp directory, and delete them at the end.

if(NOT FUZZ_TARGET OR NOT SEED_DIR OR NOT OUT_DIR)
    message(FATAL_ERROR "Usage: cmake -DFUZZ_TARGET=<target> -DSEED_DIR=<dir> -DOUT_DIR=<dir> -P RunFuzzTest.cmake")
endif()

# 1. Clean and recreate the output directory
file(REMOVE_RECURSE "${OUT_DIR}")
file(MAKE_DIRECTORY "${OUT_DIR}")

# 2. Run the fuzzer smoke test (max 10 seconds, timeout 5 seconds per input)
execute_process(
    COMMAND "${FUZZ_TARGET}" -max_total_time=10 -timeout=5 "${OUT_DIR}" "${SEED_DIR}"
    RESULT_VARIABLE res
)

# 3. Clean up the output directory and all fuzzer-generated files
file(REMOVE_RECURSE "${OUT_DIR}")

# 4. Exit with the fuzzer's result code
if(NOT res EQUAL 0)
    message(FATAL_ERROR "Fuzz test ${FUZZ_TARGET} failed with exit code ${res}")
endif()
