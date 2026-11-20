#!/bin/bash
set -eu
base="$(pwd)/.tools"
aliasbase=/home/data/t070721/tmpdir/hd_treg_v6_proseg_toolchain
if [[ ! -e "$aliasbase" ]]; then ln -s "$base" "$aliasbase"; fi
[[ "$(readlink -f "$aliasbase")" = "$base" ]]
base="$aliasbase"
mkdir -p "$base/src" "$base/rust"
for k in rustc cargo rust-std; do
 f="$k-1.88.0-x86_64-unknown-linux-gnu"
 [[ -s "$base/src/$f.tar.xz" ]] || curl -fLsS --connect-timeout 8 --max-time 120 --retry 1 "https://static.rust-lang.org/dist/2025-06-26/$f.tar.xz" -o "$base/src/$f.tar.xz"
 expected=$(python -c "import tomllib,sys; print(tomllib.load(open('logs/rust_channel.toml','rb'))['pkg'][sys.argv[1]]['target']['x86_64-unknown-linux-gnu']['xz_hash'])" "$k")
 echo "$expected  $base/src/$f.tar.xz" | sha256sum -c -
 tar -xf "$base/src/$f.tar.xz" -C "$base/src"
 "$base/src/$f/install.sh" --prefix="$base/rust" --disable-ldconfig
done
export PATH="$base/rust/bin:$PATH" CARGO_BUILD_JOBS=4 CARGO_HOME="$base/cargo"
cargo install --locked proseg --version 3.2.0 --root "$base/proseg"
"$base/proseg/bin/proseg" --version
