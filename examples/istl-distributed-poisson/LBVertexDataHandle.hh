#pragma once

#include <dune/grid/common/datahandleif.hh>

// { lb_data_handle_begin }
template <class Grid, class AssociativeContainer>
struct LBVertexDataHandle : public Dune::CommDataHandleIF<
                                LBVertexDataHandle<Grid, AssociativeContainer>,
                                typename AssociativeContainer::mapped_type> {
  LBVertexDataHandle(const std::shared_ptr<Grid> &grid,
                     AssociativeContainer &dataContainer)
      : idSet_(grid->localIdSet()), dataContainer_(dataContainer)
  {
  }

  bool contains(int dim, int codim) const
  {
    assert(dim == Grid::dimension);
    return (codim == dim); // Ony vertices have data
  }

  bool fixedSize(int dim, int codim) const
  {
    return true; // All vertices carry the same number of data items
  }

  template <class Entity> std::size_t size(const Entity &entity) const
  {
    return 1; // One data item per vertex
  }

  template <class MessageBuffer, class Entity>
  void gather(MessageBuffer &buffer, const Entity &entity) const
  {
    auto id = idSet_.id(entity);
    buffer.write(dataContainer_[id]);
  }

  template <class MessageBuffer, class Entity>
  void scatter(MessageBuffer &buffer, const Entity &entity, std::size_t n)
  {
    assert(n ==
           1); // This data handle implementations transfer only one data item.
    auto id = idSet_.id(entity);
    buffer.read(dataContainer_[id]);
  }

private:
  const typename Grid::LocalIdSet &idSet_;
  AssociativeContainer &dataContainer_;
};
// { lb_data_handle_end }