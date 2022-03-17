// { preconditioner_begin }
template <class GridView, class Matrix, class Vector>
class JacobiPreconditioner : public Dune::Preconditioner<Vector, Vector> {
public:
  // Constructor
  JacobiPreconditioner(const GridView &gridView, const Matrix &matrix)
      : gridView_(gridView), matrix_(matrix)
  {
    Vector diagonal(matrix_.N());
    for (std::size_t i = 0; i < diagonal.size(); ++i) {
      diagonal[i] = matrix_[i][i];
    }
    consistentDiagonal_ = diagonal;
    VertexDataUpdate<GridView, Vector> matrixDataHandle(gridView, diagonal,
                                                        consistentDiagonal_);

    gridView_.communicate(matrixDataHandle, All_All_Interface,
                          ForwardCommunication);
  }

  // Prepare the preconditioner
  virtual void pre(Vector &x, Vector &b) override {}

  // Apply the preconditioner
  virtual void apply(Vector &v, const Vector &r) override
  {
    auto rConsistent = r;

    VertexDataUpdate<GridView, Vector> vertexUpdateHandle(gridView_, r,
                                                          rConsistent);

    gridView_.communicate(vertexUpdateHandle,
                          InteriorBorder_InteriorBorder_Interface,
                          ForwardCommunication);

    for (std::size_t i = 0; i < matrix.N(); i++) {
      v[i] = rConsistent[i] / consistentDiagonal_[i];
    }
  }
  // Clean up
  virtual void post(Vector &x) override {}

  // Category of the preconditioner
  virtual SolverCategory::Category category() const override
  {
    return SolverCategory::sequential;
  }

private:
  const GridView gridView_;
  const Matrix &matrix_;
  Vector consistentDiagonal_;
};
// { preconditioner_end }