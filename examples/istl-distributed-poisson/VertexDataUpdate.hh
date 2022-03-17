// A DataHandle class to communicate and add vertex data
// { comm_data_handle_begin }
template <class GridView, class Vector>
struct VertexDataUpdate
    : public Dune::CommDataHandleIF<VertexDataUpdate<GridView, Vector>,
                                    typename Vector::value_type> {
  // Constructor
  VertexDataUpdate(const GridView &gridView, const Vector &userDataSend,
                   Vector &userDataReceive)
      : gridView_(gridView), userDataSend_(userDataSend),
        userDataReceive_(userDataReceive)
  {
  }

  // True if data for this codim should be communicated
  bool contains(int dim, int codim) const
  {
    return (codim == dim); // Only vertices have data
  }

  // True if data size per entity of given codim is constant
  bool fixedSize(int dim, int codim) const
  {
    return true; // All vertices carry the same number of data items
  }

  // How many objects of type DataType have to be sent for a given entity
  template <class Entiy> std::size_t size(const Entity &e) const
  {
    return 1; // One data item per vertex
  }

  // Pack user data into message buffer
  template <class MessageBuffer, class Entity>
  void gather(MessageBuffer &buffer, const Entity &entity) const
  {
    auto index = gridView_.indexSet().index(entity);
    buffer.write(userDataSend_[index]);
  }

  // Unpack user data from message buffer
  template <class MessageBuffer, class Entity>
  void scatter(MessageBuffer &buffer, const Entity &entity, std::size_t n)
  {
    assert(n == 1);
    DataType x;
    buffer.read(x);

    userDataReceive_[gridView_.indexSet().index(entity)] += x;
  }

private:
  const GridView gridView_;
  const Vector &userDataSend_;
  Vector &userDataReceive_;
};
// { comm_data_handle_end }