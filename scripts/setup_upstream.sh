#!/usr/bin/env bash
# Fetch and build the authors' OmniSketch implementation, unmodified, plus our I/O harness.
#   ext/upstream   github.com/wiegerpunter/omnisketch at commit 8dc67b5 (2023-09-04)
#   ext/jdk        Eclipse Temurin JDK 21 (only downloaded if no `javac` >= 17 is on PATH)
#   ext/build      class files
# The dependency jars shipped in the upstream repository (out/artifacts/DSCM_jar) are used.
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT=$PWD
UPSTREAM_COMMIT=8dc67b5
mkdir -p ext
if [ ! -d ext/upstream/.git ]; then
  git clone --quiet https://github.com/wiegerpunter/omnisketch "$ROOT/ext/upstream"
fi
git -C "$ROOT/ext/upstream" checkout --quiet "$UPSTREAM_COMMIT"

local_javac=$( (find ext/jdk -name javac -type f 2>/dev/null || true) | head -1)   # ext/jdk may not exist
if [ -n "$local_javac" ]; then
  JAVA_HOME_LOCAL=$ROOT/$(dirname "$(dirname "$local_javac")")
elif javac -version >/dev/null 2>&1; then   # (on macOS /usr/bin/javac may be a stub that fails)
  JAVA_HOME_LOCAL=$(dirname "$(dirname "$(readlink -f "$(command -v javac)")")")
  [ -n "${JAVA_HOME:-}" ] && [ -x "$JAVA_HOME/bin/javac" ] && JAVA_HOME_LOCAL=$JAVA_HOME
else
  case "$(uname -s)" in Darwin) os=mac ;; *) os=linux ;; esac
  case "$(uname -m)" in arm64|aarch64) arch=aarch64 ;; *) arch=x64 ;; esac
  mkdir -p ext/jdk-dl ext/jdk
  curl -sSL -o ext/jdk-dl/jdk.tar.gz \
    "https://api.adoptium.net/v3/binary/latest/21/ga/$os/$arch/jdk/hotspot/normal/eclipse"
  tar -xzf ext/jdk-dl/jdk.tar.gz -C ext/jdk
  rm -rf ext/jdk-dl
  JAVA_HOME_LOCAL=$ROOT/$(dirname "$(dirname "$(find ext/jdk -name javac -type f | head -1)")")
fi
echo "$JAVA_HOME_LOCAL" > ext/java_home.txt

CP=$(ls ext/upstream/out/artifacts/DSCM_jar/*.jar | grep -v '/DSCM\.jar$' | tr '\n' ':')
mkdir -p ext/build/upstream ext/build/harness
find ext/upstream/src/main/java -name '*.java' > ext/build/sources.txt
"$JAVA_HOME_LOCAL/bin/javac" -nowarn -d ext/build/upstream -cp "$CP" @ext/build/sources.txt
cp ext/upstream/src/main/resources/*.properties ext/build/upstream/
"$JAVA_HOME_LOCAL/bin/javac" -nowarn -d ext/build/harness -cp "ext/build/upstream:$CP" repro/omni_java/Harness.java
echo "ext/build/harness:ext/build/upstream:$CP" > ext/build/cp.txt
echo "upstream $(git -C "$ROOT/ext/upstream" rev-parse --short HEAD) built with $("$JAVA_HOME_LOCAL/bin/java" -version 2>&1 | head -1)"
