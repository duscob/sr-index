//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 6/19/2023.
//

#ifndef SRI_CONFIG_H_
#define SRI_CONFIG_H_

#include <filesystem>
#include <utility>

#include <sdsl/config.hpp>
#include <sdsl/construct_config.hpp>

#include <nlohmann/json.hpp>

namespace sri {

enum SAAlgo {
  SDSL_LIBDIVSUFSORT,
  SDSL_SE_SAIS,
  BIG_BWT
};

using JSON = nlohmann::json;

namespace conf {
constexpr std::string_view kAlphabet = "alphabet";
constexpr std::string_view kBWT = "bwt";
constexpr std::string_view kPsi = "psi";
constexpr std::string_view kBase = "base";
constexpr std::string_view kHead = "head";
constexpr std::string_view kTail = "tail";
constexpr std::string_view kPos = "pos";
constexpr std::string_view kTextPos = "textPos";
constexpr std::string_view kTextPosAsc = "textPosAsc";
constexpr std::string_view kIdx = "idx";
constexpr std::string_view kLink = "link";
}

template<uint8_t t_width>
auto createDefaultKeys() {
  using namespace conf;
  JSON keys = {
      {kAlphabet, "alphabet"},
      {kBWT, {
          {kBase, sdsl::key_bwt_trait<t_width>::KEY_BWT},
          {kHead, {
              {kPos, "bwt_run_first"},
              {kTextPos, "bwt_run_first_text_pos"},
          }},
          {kTail, {
              {kPos, "bwt_run_last"},
              {kTextPos, "bwt_run_last_text_pos"},
          }},
      },},
      {kPsi, {
          {kBase, sdsl::conf::KEY_PSI},
          {kHead, {
              {kPos, "psi_run_first"},
              {kTextPos, "psi_run_first_text_pos"},
          }},
          {kTail, {
              {kPos, "psi_run_last"},
              {kTextPos, "psi_run_last_text_pos"},
              {kTextPosAsc, {
                  {kIdx, "psi_run_last_text_pos_asc_idx"},
                  {kLink, "psi_run_last_text_pos_asc_link"},
              }},
          }},
      },},
  };

  return keys;
}

struct Config : public sdsl::cache_config {
  std::filesystem::path data_path;
  SAAlgo sa_algo = SDSL_LIBDIVSUFSORT;
  JSON keys;

  Config() = default;

  Config(const std::filesystem::path &t_data_path,
         const std::filesystem::path &t_output_dir,
         SAAlgo t_sa_algo,
         JSON t_keys = createDefaultKeys<8>())
      : data_path(t_data_path),
        sa_algo(t_sa_algo),
        keys(std::move(t_keys)),
        cache_config(false, t_output_dir, t_data_path.filename()) {
  }
};

SAAlgo toSAAlgo(const std::string &t_str) {
  static const std::map<std::string, SAAlgo> name_to_enum = {
      {"SDSL_LIBDIVSUFSORT", SDSL_LIBDIVSUFSORT},
      {"SDSL_SE_SAIS", SDSL_SE_SAIS},
      {"BIG_BWT", BIG_BWT}
  };

  return name_to_enum.at(t_str);
}

} // namespace sri

#endif //SRI_CONFIG_H_
