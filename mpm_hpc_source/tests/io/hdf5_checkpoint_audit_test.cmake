if(NOT DEFINED FIXTURE OR NOT DEFINED AUDITOR OR NOT DEFINED MODE OR
   NOT DEFINED EXPECT_SUCCESS)
  message(FATAL_ERROR "fixture, auditor, mode, and expectation are required")
endif()

set(fixture_mode "${MODE}")
if(MODE STREQUAL "legacy_allowed")
  set(fixture_mode legacy)
endif()
set(output "${CMAKE_CURRENT_BINARY_DIR}/hdf5_checkpoint_audit_${MODE}.h5")
execute_process(
  COMMAND "${FIXTURE}" "${output}" "${fixture_mode}"
  RESULT_VARIABLE fixture_status
  ERROR_VARIABLE fixture_error)
if(NOT fixture_status EQUAL 0)
  message(FATAL_ERROR "cannot create ${MODE} fixture: ${fixture_error}")
endif()

set(extra_arguments)
if(MODE STREQUAL "legacy_allowed")
  list(APPEND extra_arguments --allow-legacy)
endif()
execute_process(
  COMMAND "${AUDITOR}" --hdf5 "${output}" --jsonl ${extra_arguments}
  RESULT_VARIABLE audit_status
  OUTPUT_VARIABLE audit_output
  ERROR_VARIABLE audit_error)
file(REMOVE "${output}")

if(EXPECT_SUCCESS)
  if(NOT audit_status EQUAL 0)
    message(FATAL_ERROR "audit unexpectedly failed: ${audit_error}")
  endif()
  string(REGEX MATCH "\"schema\":\"mpm-hdf5-checkpoint-audit-v1\""
         schema_match "${audit_output}")
  if(NOT schema_match)
    message(FATAL_ERROR "audit omitted schema: ${audit_output}")
  endif()
  string(REGEX MATCHALL "\"record_type\":\"particle\""
         particle_rows "${audit_output}")
  list(LENGTH particle_rows particle_count)
  if(NOT particle_count EQUAL 3)
    message(FATAL_ERROR "expected 3 JSONL particles, got ${particle_count}")
  endif()
else()
  if(audit_status EQUAL 0)
    message(FATAL_ERROR "audit unexpectedly passed: ${audit_output}")
  endif()
endif()
