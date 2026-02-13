# frozen_string_literal: true

require "gelfgren"

RSpec.configure do |config|
  # Enable flags like --only-failures and --next-failure
  config.example_status_persistence_file_path = ".rspec_status"

  # Disable RSpec exposing methods globally on `Module` and `main`
  config.disable_monkey_patching!

  config.expect_with :rspec do |c|
    c.syntax = :expect
  end

  # Custom matchers
  RSpec::Matchers.define :be_close_to do |expected, tolerance: 1e-10|
    match do |actual|
      (actual - expected).abs < tolerance
    end

    failure_message do |actual|
      "expected #{actual} to be within #{tolerance} of #{expected}, " \
        "but difference was #{(actual - expected).abs}"
    end
  end
end
