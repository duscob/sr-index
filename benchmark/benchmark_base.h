//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 12/25/25.
//

#pragma once

#include <functional>

#include "sr-index/alphabet.h"
#include "sr-index/index_base.h"


#ifndef SRI_BENCHMARK_DATA_WIDTH
#define SRI_BENCHMARK_DATA_WIDTH 8
#endif


using ExternalGenericStorage = std::reference_wrapper<sri::GenericStorage>;

using Alphabet = sri::Alphabet<SRI_BENCHMARK_DATA_WIDTH>;
