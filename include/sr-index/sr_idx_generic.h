#pragma once

namespace sri {

template <typename TSrIdx, uint32_t TSubsampleRate>
class SrIdxGeneric : public TSrIdx {
 public:
  using Base = TSrIdx;

  SrIdxGeneric(const typename TSrIdx::Storage& t_storage) : TSrIdx(t_storage, TSubsampleRate) {};
  SrIdxGeneric() : TSrIdx(TSubsampleRate) {};
};

//~~~~~~~

template <typename TSrIdx, uint32_t TSubsampleRate>
void construct(SrIdxGeneric<TSrIdx, TSubsampleRate>& t_index, const std::string& t_data_path, Config& t_config) {
  using Index = SrIdxGeneric<TSrIdx, TSubsampleRate>;
  construct(dynamic_cast<typename Index::Base&>(t_index), t_data_path, t_config);

  t_index.load(t_config);
}

}  // namespace sri