#define CATCH_CONFIG_NO_POSIX_SIGNALS
#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <stdexcept>

#include "solvers/prescribed_pressure_validation.h"

TEST_CASE("Pressure database writer rejects a history beyond the analysis") {
  using mpm::prescribed_pressure::detail::validate_write_max_step;

  REQUIRE_NOTHROW(validate_write_max_step(0, 0));
  REQUIRE_NOTHROW(validate_write_max_step(130, 104000));
  REQUIRE_NOTHROW(validate_write_max_step(104000, 104000));
  REQUIRE_THROWS_AS(validate_write_max_step(104001, 104000),
                    std::invalid_argument);
}
