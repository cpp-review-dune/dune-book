// Scalar product for pairs of additive and consistent vectors
template <class GridView, class Vector>
// { scalar_product_begin }
class AdditiveScalarProduct : public ScalarProduct<Vector> {
public:
  // Constructor
  AdditiveScalarProduct(const GridView &gridView) : gridView_(gridView) {}

  // Dot product of a consistent vector x and an additive vector y,
  // or viceversa
  virtual field_type dot(const Vector &x, const Vector &y) const override
  {
    return gridView_.comm().sum(x.dot(y));
  }

  // Norm of a vector given in additive representation
  virtual real_type norm(const Vector &x) const override
  {
    // Construct consistent representation of x
    auto xConsistent = x;

    VertexDataUpdate<GridView, Vector> vertexUpdateHandle(gridView_, x,
                                                          xConsistent);

    gridView_.communicate(vertexUpdateHandle,
                          InteriorBorder_InteriorBorder_Interface,
                          ForwardCommunication);

    // Local scalar product of x with itself
    auto localNorm2 = x.dot(xConsistent);

    // Sum over all subdomains
    return std::sqrt(gridView_.comm().sum(localNorm2));
  }

  // Scalar product, linear operator, and preconditioner must be
  // of the same category
  virtual SolverCategory::Category category() const override
  {
    return SolverCategory::sequential;
  }

private:
  const GridView gridView_;
};
// { scalar_product_end }