# frozen_string_literal: true

Gem::Specification.new do |spec|
  spec.name = "gelfgren"
  spec.version = "0.1.0"
  spec.authors = ["Gelfgren Contributors"]
  spec.email = ["your.email@example.com"]

  spec.summary = "High-performance numerical computing library with piecewise rational interpolation"
  spec.description = <<-DESC
    Gelfgren is a numerical computing library implementing piecewise rational interpolation
    methods based on Jan Gelfgren's 1975 research. It provides Bernstein polynomials,
    rational functions, Padé approximants, and boundary value problem solvers.
    Written in Rust with Ruby bindings via Magnus.
  DESC
  spec.homepage = "https://github.com/yourusername/gelfgren"
  spec.license = "MIT OR Apache-2.0"
  spec.required_ruby_version = ">= 3.0.0"

  spec.metadata["homepage_uri"] = spec.homepage
  spec.metadata["source_code_uri"] = "https://github.com/yourusername/gelfgren"
  spec.metadata["changelog_uri"] = "https://github.com/yourusername/gelfgren/blob/main/CHANGELOG.md"
  spec.metadata["documentation_uri"] = "https://yourusername.github.io/gelfgren"

  # Specify which files should be added to the gem
  spec.files = Dir[
    "lib/**/*.rb",
    "ext/**/*.{rs,toml,rb}",
    "README.md",
    "LICENSE*"
  ]

  spec.bindir = "exe"
  spec.executables = spec.files.grep(%r{\Aexe/}) { |f| File.basename(f) }
  spec.require_paths = ["lib"]
  spec.extensions = ["ext/gelfgren/extconf.rb"]

  # Runtime dependencies
  spec.add_dependency "rb_sys", "~> 0.9"

  # Development dependencies
  spec.add_development_dependency "rake", "~> 13.0"
  spec.add_development_dependency "rake-compiler", "~> 1.2"
  spec.add_development_dependency "rspec", "~> 3.12"
  spec.add_development_dependency "rubocop", "~> 1.50"

  # For building the Rust extension
  spec.metadata["cargo_toml"] = "ext/gelfgren/Cargo.toml"
end
