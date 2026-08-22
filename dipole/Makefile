# =====================================================
# Compilador e flags
# =====================================================

CXX := g++

# -MMD -MP gera .d com as dependencias de header. Sem isso, editar
# include/*.hpp NAO dispara recompilacao -- e a campanha edita headers.
CXXFLAGS := -O3 -std=c++17 -Wall -Wextra -fopenmp -MMD -MP -Iinclude -I$(CONDA_PREFIX)/include $(shell lhapdf-config --cflags)
LDFLAGS  := -fopenmp -L$(CONDA_PREFIX)/lib $(shell lhapdf-config --libs)


# =====================================================
# Diretórios
# =====================================================

SRC_DIR   := src
BUILD_DIR := build
DATA_DIR  := data
PLOT_DIR  := plots


# =====================================================
# Parâmetros globais
# =====================================================

# Energia do neutrino em GeV.
# Pode ser sobrescrita na chamada:
#
# make plot-tau-profile-collapsar ENU=1e12
#
ENU := 1e12

# Faixa padrão para calcular sigma_nuN(E_nu)
LOGEMIN := 3
LOGEMAX := 14
NE_SIGMA := 300

# =====================================================
# Modelo de dipolo
#
# Opções:
# GBW
# IIM
#
# Pode sobrescrever:
#
# make sigma-nuN DIPOLE_MODEL=IIM
# =====================================================

DIPOLE_MODEL := GBW

# =====================================================
# F3

USE_F3 := 1
# PDF_SET := NNPDF31_nlo_as_0118
PDF_SET := CT10nlo
BEAM := nu

# =====================================================
# Perfil de Collapsar usado nos cálculos
# =====================================================

# PROFILE_NAME ?= benchmark_01_radiating_sphere
PROFILE_NAME ?= benchmark_02_collapsar_disk_funnel
#PROFILE_NAME ?= benchmark_03_uhe_attenuation

PROFILE_FILE := data/$(PROFILE_NAME).dat

# =====================================================
# Executáveis
# =====================================================

TARGET             := $(BUILD_DIR)/dipole
SCAN_X             := $(BUILD_DIR)/scan_x
SIGMA_NUN          := $(BUILD_DIR)/sigma_nuN
PROFILE_GENERATOR := $(BUILD_DIR)/generate_profile
TAU                := $(BUILD_DIR)/optical_depth
TAU_PROFILE        := $(BUILD_DIR)/tau_profile
COLLAPSAR_PROFILE  := $(BUILD_DIR)/generate_collapsar_profile
SCAN_NEUTRINOSPHERE := $(BUILD_DIR)/scan_neutrinosphere
BENCHMARK_PROFILES := $(BUILD_DIR)/generate_benchmark_profiles
VALIDATE_HERA      := $(BUILD_DIR)/validate_hera


# =====================================================
# Fontes principais do executável dipole
# =====================================================

SOURCES := \
	$(SRC_DIR)/main.cpp \
	$(SRC_DIR)/wavefunctions.cpp \
	$(SRC_DIR)/dipole_models.cpp \
	$(SRC_DIR)/integrals.cpp

OBJECTS := $(SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)


# =====================================================
# Atalhos rápidos
# =====================================================

gbw:
	$(MAKE) -B collapsar-all DIPOLE_MODEL=GBW

iim:
	$(MAKE) -B collapsar-all DIPOLE_MODEL=IIM

# =====================================================
# Alvos que não correspondem a arquivos físicos
# =====================================================

.PHONY: \
	all clean dirs run \
	wavefunctions integrand sigma \
	plots plot-wavefunctions plot-integrand plot-sigma \
	scan-x plot-scan-x \
	sigma-nuN plot-sigma-nuN \
	profile plot-profile \
	tau plot-tau \
	tau-profile plot-tau-profile \
	collapsar-profile tau-profile-collapsar plot-tau-profile-collapsar \
	collapsar-all \
	scan-neutrinosphere plot-neutrinosphere-scan \
	benchmark-profiles \
	gbw iim


# =====================================================
# Construção básica
# =====================================================

all: dirs $(TARGET)

dirs:
	mkdir -p $(BUILD_DIR) $(DATA_DIR) $(PLOT_DIR)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


# =====================================================
# Execução básica do código principal
# =====================================================

run: all
	$(TARGET)

wavefunctions: all
	$(TARGET) --mode wavefunctions --current NC --flavor u --Q2 10.0

integrand: all
	$(TARGET) --mode integrand --current NC --flavor u --x 1e-5 --Q2 10.0

sigma: all
	$(TARGET) --mode sigma --current NC --flavor u --x 1e-5 --Q2 10.0


# =====================================================
# Plots básicos
# =====================================================

plots: plot-wavefunctions

plot-wavefunctions: wavefunctions
	python $(PLOT_DIR)/plot_wavefunctions.py

plot-integrand: integrand
	python $(PLOT_DIR)/plot_integrand.py

plot-sigma: sigma
	python $(PLOT_DIR)/plot_sigma.py


# =====================================================
# Scan em x para F2(x,Q2)
# =====================================================

scan-x: dirs $(SCAN_X)
	LD_LIBRARY_PATH=$(CONDA_PREFIX)/lib:$$LD_LIBRARY_PATH $(SCAN_X) \
		--use-F3 $(USE_F3) \
		--pdf-set $(PDF_SET) \
		--beam $(BEAM)

$(SCAN_X): $(BUILD_DIR)/scan_x.o \
           $(BUILD_DIR)/wavefunctions.o \
           $(BUILD_DIR)/dipole_models.o \
           $(BUILD_DIR)/integrals.o \
           $(BUILD_DIR)/weak_structure_functions.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILD_DIR)/scan_x.o: $(SRC_DIR)/scan_x.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

plot-scan-x: scan-x
	python $(PLOT_DIR)/plot_scan_x_models.py


# =====================================================
# Seção de choque sigma_nuN(E_nu)
# =====================================================

sigma-nuN: dirs $(SIGMA_NUN)
	LD_LIBRARY_PATH=$(CONDA_PREFIX)/lib:$$LD_LIBRARY_PATH $(SIGMA_NUN) \
		--model $(DIPOLE_MODEL) \
		--logEmin $(LOGEMIN) \
		--logEmax $(LOGEMAX) \
		--NE $(NE_SIGMA) \
		--use-F3 $(USE_F3) \
		--pdf-set $(PDF_SET) \
		--beam $(BEAM)

$(SIGMA_NUN): $(BUILD_DIR)/sigma_nuN.o $(BUILD_DIR)/wavefunctions.o $(BUILD_DIR)/dipole_models.o $(BUILD_DIR)/integrals.o $(BUILD_DIR)/structure_table.o $(BUILD_DIR)/sigma_nuN_core.o $(BUILD_DIR)/weak_structure_functions.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILD_DIR)/sigma_nuN.o: $(SRC_DIR)/sigma_nuN.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/weak_structure_functions.o: $(SRC_DIR)/weak_structure_functions.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# =====================================================
# F1 - validacao do limite EM contra o F2 do HERA
#
# GBW e bCGC sao ajustes ao HERA. Se o limite EM nao reproduz
# esse dado, nada a jusante vale. Ver CAMPANHA_CORRECAO.md.
#
#   make validate-hera                 # grade default
#   make validate-hera NR=800 NZ=400   # teste de convergencia
# =====================================================

NR := 400
NZ := 200
HERA_FLAVORS := uds
HERA_MQ := -1

validate-hera: dirs $(VALIDATE_HERA)
	LD_LIBRARY_PATH=$(CONDA_PREFIX)/lib:$$LD_LIBRARY_PATH $(VALIDATE_HERA) \
		--Nr $(NR) --Nz $(NZ) \
		--flavors $(HERA_FLAVORS) \
		--mq $(HERA_MQ)

$(VALIDATE_HERA): $(BUILD_DIR)/validate_hera.o $(BUILD_DIR)/wavefunctions.o $(BUILD_DIR)/dipole_models.o $(BUILD_DIR)/integrals.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILD_DIR)/validate_hera.o: $(SRC_DIR)/validate_hera.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

TEST_QUAD := $(BUILD_DIR)/test_quadrature

test-quadrature: dirs $(TEST_QUAD)
	LD_LIBRARY_PATH=$(CONDA_PREFIX)/lib:$$LD_LIBRARY_PATH $(TEST_QUAD)

$(TEST_QUAD): $(BUILD_DIR)/test_quadrature.o $(BUILD_DIR)/wavefunctions.o $(BUILD_DIR)/dipole_models.o $(BUILD_DIR)/integrals.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILD_DIR)/test_quadrature.o: $(SRC_DIR)/test_quadrature.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

TEST_TABLE := $(BUILD_DIR)/test_table

test-table: dirs $(TEST_TABLE)
	LD_LIBRARY_PATH=$(CONDA_PREFIX)/lib:$$LD_LIBRARY_PATH $(TEST_TABLE)

$(TEST_TABLE): $(BUILD_DIR)/test_table.o $(BUILD_DIR)/wavefunctions.o $(BUILD_DIR)/dipole_models.o $(BUILD_DIR)/integrals.o $(BUILD_DIR)/structure_table.o $(BUILD_DIR)/sigma_nuN_core.o $(BUILD_DIR)/weak_structure_functions.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILD_DIR)/test_table.o: $(SRC_DIR)/test_table.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

plot-sigma-nuN: sigma-nuN
	python $(PLOT_DIR)/plot_sigma_nuN.py


# =====================================================
# Perfil FLASH-like simples
# =====================================================

profile: dirs $(PROFILE_GENERATOR)
	$(PROFILE_GENERATOR)

$(PROFILE_GENERATOR): $(BUILD_DIR)/generate_profile.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILD_DIR)/generate_profile.o: $(SRC_DIR)/generate_profile.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

plot-profile: profile
	python $(PLOT_DIR)/plot_profile.py


# =====================================================
# Profundidade óptica total tau(E_nu)
# =====================================================

tau: sigma-nuN benchmark-profiles dirs $(TAU)
	$(TAU) \
		--profile $(PROFILE_FILE) \
		--sigma-table data/sigma_nuN_CC_$(DIPOLE_MODEL).dat \
		--output data/optical_depth_CC_$(DIPOLE_MODEL).dat

$(TAU): $(BUILD_DIR)/optical_depth.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILD_DIR)/optical_depth.o: $(SRC_DIR)/optical_depth.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

plot-tau: tau
	python $(PLOT_DIR)/plot_optical_depth.py $(DIPOLE_MODEL) $(PROFILE_NAME)

# =====================================================
# Perfil radial de tau para perfil FLASH-like
# Usa sigma_nuN(E_nu) calculada pelo método de dipolo.
# =====================================================

tau-profile: sigma-nuN profile dirs $(TAU_PROFILE)
	$(TAU_PROFILE) \
		--profile $(PROFILE_FILE) \
		--sigma-table data/sigma_nuN_CC_$(DIPOLE_MODEL).dat \
		--output data/tau_profile.dat \
		--Enu $(ENU)

$(TAU_PROFILE): $(BUILD_DIR)/tau_profile.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILD_DIR)/tau_profile.o: $(SRC_DIR)/tau_profile.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

plot-tau-profile: tau-profile
	python $(PLOT_DIR)/plot_tau_profile.py data/tau_profile.dat $(PROFILE_NAME) $(PROFILE_NAME)


# =====================================================
# Perfil de collapsar
# =====================================================

collapsar-profile: dirs $(COLLAPSAR_PROFILE)
	$(COLLAPSAR_PROFILE)

$(COLLAPSAR_PROFILE): $(BUILD_DIR)/generate_collapsar_profile.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILD_DIR)/generate_collapsar_profile.o: $(SRC_DIR)/generate_collapsar_profile.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


# =====================================================
# Perfil radial de tau para collapsar
# Usa sigma_nuN(E_nu) calculada pelo método de dipolo.
#
# Uso:
# make plot-tau-profile-collapsar ENU=1e9
# make plot-tau-profile-collapsar ENU=1e12
# =====================================================

tau-profile-collapsar: sigma-nuN benchmark-profiles dirs $(TAU_PROFILE)
	$(TAU_PROFILE) \
		--profile $(PROFILE_FILE) \
		--sigma-table data/sigma_nuN_CC_$(DIPOLE_MODEL).dat \
		--output data/tau_profile_collapsar.dat \
		--Enu $(ENU)

plot-tau-profile-collapsar: tau-profile-collapsar
	python $(PLOT_DIR)/plot_tau_profile.py data/tau_profile_collapsar.dat collapsar


# =====================================================
# Pipeline completo de collapsar
# Gera perfil, calcula sigma_nuN, calcula tau(r), plota.
#
# Uso:
# make collapsar-all ENU=1e12
# =====================================================

collapsar-all: plot-tau-profile-collapsar


# =====================================================
# Scan sistemático da neutrinosfera UHE
# Gera:
# data/neutrinosphere_scan.dat
# data/tau_profiles_scan.dat
# =====================================================

scan-neutrinosphere: sigma-nuN benchmark-profiles dirs $(SCAN_NEUTRINOSPHERE)
	$(SCAN_NEUTRINOSPHERE) --model $(DIPOLE_MODEL) --profile $(PROFILE_FILE)

$(SCAN_NEUTRINOSPHERE): $(BUILD_DIR)/scan_neutrinosphere.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILD_DIR)/scan_neutrinosphere.o: $(SRC_DIR)/scan_neutrinosphere.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

plot-neutrinosphere-scan: scan-neutrinosphere
	python $(PLOT_DIR)/plot_neutrinosphere_scan.py $(DIPOLE_MODEL) $(PROFILE_NAME)



# =====================================================
# Limpeza
# =====================================================

clean:
	rm -rf $(BUILD_DIR)/*.o \
	       $(TARGET) \
	       $(SCAN_X) \
	       $(SIGMA_NUN) \
	       $(PROFILE_GENERATOR) \
	       $(TAU) \
	       $(TAU_PROFILE) \
	       $(COLLAPSAR_PROFILE) \
	       $(SCAN_NEUTRINOSPHERE)

# =====================================================
# Benchmark profiles para testes DIS / tau(E)
# =====================================================

$(BENCHMARK_PROFILES): $(SRC_DIR)/generate_benchmark_profiles.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

benchmark-profiles: dirs $(BENCHMARK_PROFILES)
	$(BENCHMARK_PROFILES)

# =====================================================
# Dependencias de header geradas por -MMD -MP
# =====================================================

-include $(wildcard $(BUILD_DIR)/*.d)
