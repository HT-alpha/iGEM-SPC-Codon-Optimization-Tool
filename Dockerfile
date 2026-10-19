# DNA Design Suite. The build stage compiles the Rust core. The runtime
# image is Python plus that wheel: a host with Docker does not install Rust.
# Toolchain pin matches rust-toolchain.toml.
# Ubuntu is the base so apt can install the compiler.

FROM ubuntu:24.04 AS build

ENV DEBIAN_FRONTEND=noninteractive \
    RUSTUP_HOME=/usr/local/rustup \
    CARGO_HOME=/usr/local/cargo \
    PATH=/usr/local/cargo/bin:$PATH \
    CARGO_BUILD_JOBS=2

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        pkg-config \
        python-is-python3 \
        python3 \
        python3-dev \
        python3-pip \
        python3-venv \
    && rm -rf /var/lib/apt/lists/*

RUN curl --proto '=https' --tlsv1.2 -fsSL https://sh.rustup.rs \
    | sh -s -- -y --profile minimal --default-toolchain 1.98.1 --component rustfmt,clippy

WORKDIR /src
COPY . .
RUN python tools/stage_wheel_examples.py \
    && python3 -m venv /opt/build-venv \
    && /opt/build-venv/bin/pip install --no-cache-dir "maturin>=1.5,<2" \
    && /opt/build-venv/bin/maturin build --release --interpreter python3 --out /wheels

FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PATH=/opt/dna-design-suite/venv/bin:$PATH

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        python3 \
        python3-venv \
    && rm -rf /var/lib/apt/lists/* \
    && python3 -m venv /opt/dna-design-suite/venv

COPY --from=build /wheels /tmp/wheels
RUN /opt/dna-design-suite/venv/bin/pip install --no-cache-dir /tmp/wheels/*.whl \
    && rm -rf /tmp/wheels

COPY examples/ /opt/dna-design-suite/examples/

WORKDIR /work
EXPOSE 8000
ENTRYPOINT ["dna-design-suite"]
