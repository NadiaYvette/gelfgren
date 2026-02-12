Name:           gelfgren
Version:        0.1.0
Release:        1%{?dist}
Summary:        Piecewise rational interpolation library

License:        MIT OR Apache-2.0
URL:            https://github.com/yourusername/gelfgren
Source0:        %{name}-%{version}.tar.gz

BuildRequires:  rust >= 1.70
BuildRequires:  cargo
BuildRequires:  gcc
BuildRequires:  make

%description
Gelfgren is a high-performance numerical computing library implementing
piecewise rational interpolation methods based on Jan Gelfgren's 1975 research.
Written in Rust with bindings for 17 programming languages.

%package        devel
Summary:        Development files for %{name}
Requires:       %{name}%{?_isa} = %{version}-%{release}

%description    devel
The %{name}-devel package contains libraries and header files for
developing applications that use %{name}.

%prep
%autosetup

%build
cargo build --release --manifest-path gelfgren-core/Cargo.toml
cargo build --release --manifest-path gelfgren-ffi/Cargo.toml
cd gelfgren-ffi
cbindgen --config cbindgen.toml --output ../include/gelfgren.h

%install
# Install library
mkdir -p %{buildroot}%{_libdir}
install -m 755 target/release/libgelfgren.so %{buildroot}%{_libdir}/libgelfgren.so.%{version}
ln -s libgelfgren.so.%{version} %{buildroot}%{_libdir}/libgelfgren.so.0
ln -s libgelfgren.so.0 %{buildroot}%{_libdir}/libgelfgren.so

# Install headers
mkdir -p %{buildroot}%{_includedir}
install -m 644 include/gelfgren.h %{buildroot}%{_includedir}/

# Install pkgconfig
mkdir -p %{buildroot}%{_libdir}/pkgconfig
cat > %{buildroot}%{_libdir}/pkgconfig/gelfgren.pc << EOF
prefix=%{_prefix}
exec_prefix=%{_exec_prefix}
libdir=%{_libdir}
includedir=%{_includedir}

Name: Gelfgren
Description: Piecewise rational interpolation library
Version: %{version}
Libs: -L\${libdir} -lgelfgren
Cflags: -I\${includedir}
EOF

%ldconfig_scriptlets

%files
%license LICENSE-MIT LICENSE-APACHE
%doc README.md
%{_libdir}/libgelfgren.so.0*
%{_libdir}/libgelfgren.so.%{version}

%files devel
%{_includedir}/gelfgren.h
%{_libdir}/libgelfgren.so
%{_libdir}/pkgconfig/gelfgren.pc

%changelog
* Wed Feb 13 2026 Nadia Chambers <nadia.chambers@iohk.io> - 0.1.0-1
- Initial RPM package
- Support for Bernstein polynomials, rational functions, and Padé approximants
- Bindings for C, C++, Python, Java, R, Ruby, Fortran, Haskell, Mercury,
  OCaml, Julia, Go, Standard ML, Common Lisp, Scheme, and Prolog
