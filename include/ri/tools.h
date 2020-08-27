//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/27/20.
//

#ifndef RI_TOOLS_H_
#define RI_TOOLS_H_

#include <cstddef>

namespace ri {

template<typename TContainer>
class RandomAccessForContainer {
 public:
  RandomAccessForContainer(const TContainer &t_container) : container_{t_container} {}

  auto operator()(std::size_t i) const {
    return container_.get()[i];
  }

 private:
  TContainer container_;
};

template<typename TContainer>
auto buildRandomAccessForContainer(const TContainer &t_container) {
  return RandomAccessForContainer<TContainer>(t_container);
}

template<typename TContainer1, typename TContainer2>
class RandomAccessForTwoContainers {
 public:
  RandomAccessForTwoContainers(const TContainer1 &t_container1, const TContainer2 &t_container2)
      : container1_{t_container1}, container2_{t_container2} {
  }

  auto operator()(std::size_t i) const{
    return std::make_pair(container1_.get()[i], container2_.get()[i]);
  }

 private:
  TContainer1 container1_;
  TContainer2 container2_;
};

template<typename TContainer1, typename TContainer2>
auto buildRandomAccessForTwoContainers(const TContainer1 &t_container1, const TContainer2 &t_container_2) {
  return RandomAccessForTwoContainers<TContainer1, TContainer2>(t_container1, t_container_2);
}

}

#endif //RI_TOOLS_H_
