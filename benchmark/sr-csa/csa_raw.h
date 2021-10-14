//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 10/12/21.
//

#ifndef SRI_BENCHMARK_SR_CSA_CSA_RAW_H_
#define SRI_BENCHMARK_SR_CSA_CSA_RAW_H_

template<uint8_t t_width = 8,
    typename TAlphabet = sdsl::byte_alphabet,
    typename TPsiRLE = sri::PsiCoreRLE<>,
    typename TSA = sdsl::int_vector<>>
class CSARaw : public CSA<t_width, TAlphabet, TPsiRLE> {
 public:

  using BaseClass = CSA<t_width, TAlphabet, TPsiRLE>;

  explicit CSARaw(std::reference_wrapper<ExternalStorage> t_storage) : BaseClass(t_storage) {}

  using typename BaseClass::size_type;

  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v, const std::string &name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

    size_type written_bytes = 0;
    written_bytes += this->template serializeItem<TAlphabet>(
        sri::key_trait<t_width>::KEY_ALPHABET, out, child, "alphabet");

    written_bytes += this->template serializeItem<TPsiRLE>(sdsl::conf::KEY_PSI, out, child, "psi");

    written_bytes += this->template serializeItem<TSA>(sdsl::conf::KEY_SA, out, child, "sa");

    return written_bytes;
  }

 protected:

  using typename BaseClass::TSource;
  using typename BaseClass::TFnComputeSAValues;
  TFnComputeSAValues constructComputeSAValues(TSource &t_source) override {
    auto cref_sa = this->template loadItem<TSA>(sdsl::conf::KEY_SA, t_source);

    auto compute_sa_values = [cref_sa](const auto &tt_range, const auto &, auto tt_report) {
      auto[first, last] = tt_range;

      while (first <= last) {
        tt_report(cref_sa.get()[first++]);
      }
    };

    return compute_sa_values;
  }

};

#endif //SRI_BENCHMARK_SR_CSA_CSA_RAW_H_
