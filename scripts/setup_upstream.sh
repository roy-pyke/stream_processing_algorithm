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

if [ -x ext/jdk/bin/javac ]; then
  JAVA_HOME_LOCAL=$ROOT/ext/jdk
elif command -v javac >/dev/null 2>&1; then
  JAVA_HOME_LOCAL=$(dirname "$(dirname "$(readlink -f "$(command -v javac)")")")
else
  mkdir -p ext/jdk-dl ext/jdk
  curl -sSL -o ext/jdk-dl/jdk.tar.gz \
    "https://api.adoptium.net/v3/binary/latest/21/ga/linux/x64/jdk/hotspot/normal/eclipse"
  tar -xzf ext/jdk-dl/jdk.tar.gz -C ext/jdk --strip-components=1
  rm -rf ext/jdk-dl
  JAVA_HOME_LOCAL=$ROOT/ext/jdk
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
