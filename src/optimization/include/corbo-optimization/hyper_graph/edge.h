/*********************************************************************
 *
 *  Software License Agreement
 *
 *  Copyright (c) 2020,
 *  TU Dortmund - Institute of Control Theory and Systems Engineering.
 *  All rights reserved.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  Authors: Christoph Rösmann
 *********************************************************************/

#ifndef SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_EDGE_H_
#define SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_EDGE_H_

#include <corbo-core/types.h>
#include <corbo-optimization/hyper_graph/edge_interface.h>
#include <corbo-optimization/hyper_graph/vector_vertex.h>

#include <corbo-core/utilities.h>

#include <array>
#include <memory>

namespace corbo {

/**
 * @brief Templated base edge class that stores an arbitary number of value
 *
 * @ingroup optimization hyper-graph
 *
 * The value dimension is specified statically via a template parameter.
 * The dimension has to be known at compile time and is set by the template parameter
 * \c D.
 *
 * Connected vertex types are defined with a number of class templates (Vertices...).
 * Internally, those vertices are stored
 * according to their order of appearence in a vertex container _vertices.
 *
 * Consequently, in the computeValues() method, vertices can be accessed by
 * _vertices[idx] and statically casted to their particular type if necessary.
 *
 * @warning Vertices must remain valid as long as the Edge is in use!
 *
 * @see EdgeInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 *
 * @tparam D Dimension of the edge (value vector).
 * @tparam Vertices... An arbitary number of vertex types that are attached to this edge
 */
// template <class... Vertices>
// class Edge : public BaseEdge
//{
// public:
//    using Ptr  = std::shared_ptr<Edge>;
//    using UPtr = std::unique_ptr<Edge>;

//    //! Return number of vertices at compile-time
//    static constexpr const int numVerticesCompileTime = sizeof...(Vertices);
//    //! Return edge dimension at compile-time
//    static constexpr const int dimensionCompileTime = D;

//    //! Typedef to represent the vertex container.
//    using VertexContainer = std::array<VertexInterface*, numVerticesCompileTime>;

//    // delete default constructor
//    Edge() = delete;

//    /**
//     * @brief Construct edge by providing connected vertices
//     *
//     * The order must match the order of template arguments provided to the class.
//     * @warning Vertices must remain valid as long as the Edge is in use!
//     */
//    template <class... VerticesT>
//    explicit Edge(VerticesT&... args) : BaseEdge(), _vertices({&args...})
//    {
//        // Check if types and number of template parameters match class template paremters for vertices
//        static_assert(util::variadic_temp_equal<std::tuple<Vertices...>, std::tuple<VerticesT...>>::value,
//                      "BaseEdge(): Number and types of vertices passed via the constructor does not match number of class template parameters");
//    }  //!< Construct edge by passing all vertices references.

//    // implements interface method
//    int getDimension() const override { return D; }
//    // implements interface method
//    int getNumVertices() const override { return numVerticesCompileTime; }

//    // implements interface method
//    int verticesDimension() const override
//    {
//        int n = 0;
//        for (VertexInterface* vertex : _vertices) n += vertex->getDimension();
//        return n;
//    }

//    // implements interface method
//    bool isLinear() const override = 0;
//    // implements interface method
//    bool providesJacobian() const override { return false; }

//    // implements interface method
//    void computeValues(double* values) override = 0;

//    // implements interface method
//    VertexInterface* getVertexRaw(int idx) override
//    {
//        assert(idx < getNumVertices());
//        return _vertices[idx];
//    }

// protected:
//    const VertexContainer _vertices;  //!< Vertex container

// public:
//    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
//};

/**
 * @brief Templated base edge class that stores an arbitary number of value
 *
 * @ingroup optimization hyper-graph
 *
 * Template specialization for edges with a dynamic value dimension.
 * Refer to the description of the non-specialized class.
 *
 * @warning Vertices must remain valid as long as the Edge is in use!
 *
 * @see EdgeInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 *
 * @tparam D Dimension of the edge (value vector).
 */
template <class... Vertices>
class Edge : public BaseEdge
{
 public:
    using Ptr      = std::shared_ptr<Edge>;
    using ConstPtr = std::shared_ptr<const Edge>;
    using UPtr     = std::unique_ptr<Edge>;

    //! Return number of vertices at compile-time
    static constexpr const int numVerticesCompileTime = sizeof...(Vertices);

    //! Typedef to represent the vertex container.
    using VertexContainer = std::array<VertexInterface*, numVerticesCompileTime>;

    // delete default constructor
    Edge() = delete;

    /**
     * @brief Construct edge by providing connected vertices
     *
     * The order must match the order of template arguments provided to the class.
     * @warning Vertices must remain valid as long as the Edge is in use!
     */
    template <class... VerticesT>
    explicit Edge(VerticesT&... args) : BaseEdge(), _vertices({&args...})
    {
        // Check if types and number of template parameters match class template paremters for vertices
        static_assert(util::variadic_temp_equal<std::tuple<Vertices...>, std::tuple<VerticesT...>>::value,
                      "Edge(): Number and types of vertices passed via the constructor does not match number of class template parameters");
    }  //!< Construct edge by passing all vertices references.

    // implements interface method
    int getDimension() const override = 0;
    // implements interface method
    int getNumVertices() const override { return numVerticesCompileTime; }

    // implements interface method
    int verticesDimension() const override
    {
        int n = 0;
        for (VertexInterface* vertex : _vertices) n += vertex->getDimension();
        return n;
    }

    // implements interface method
    bool isLinear() const override = 0;
    // implements interface method
    bool providesJacobian() const override { return false; }

    // implements interface method
    void computeValues(Eigen::Ref<Eigen::VectorXd> values) override = 0;

    // implements interface method
    VertexInterface* getVertexRaw(int idx) override
    {
        assert(idx < getNumVertices());
        return _vertices[idx];
    }

    // implements interface method
    const VertexInterface* getVertex(int idx) const override
    {
        assert(idx < getNumVertices());
        return _vertices[idx];
    }

 protected:
    const VertexContainer _vertices;  //!< Vertex container

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template <>
class Edge<> : public BaseEdge
{
 public:
    using Ptr      = std::shared_ptr<Edge>;
    using ConstPtr = std::shared_ptr<const Edge>;
    using UPtr     = std::unique_ptr<Edge>;

    //! Typedef to represent the vertex container.
    using VertexContainer = std::vector<VertexInterface*>;

    // delete default constructor
    Edge() = delete;

    explicit Edge(int num_vertices) { resizeVertexContainer(num_vertices); }

    //! Set number \c n of vertices attached to this edge.
    void resizeVertexContainer(int n) { _vertices.resize(n); }

    void setVertex(int idx, VertexInterface& vertex)
    {
        assert(idx < (int)_vertices.size());
        _vertices[idx] = &vertex;
    }

    // implements interface method
    int getDimension() const override = 0;
    // implements interface method
    int getNumVertices() const override { return _vertices.size(); }

    // implements interface method
    int verticesDimension() const override
    {
        int n = 0;
        for (VertexInterface* vertex : _vertices) n += vertex->getDimension();
        return n;
    }

    // implements interface method
    bool isLinear() const override = 0;
    // implements interface method
    bool providesJacobian() const override { return false; }

    // implements interface method
    void computeValues(Eigen::Ref<Eigen::VectorXd> values) override = 0;

    // implements interface method
    VertexInterface* getVertexRaw(int idx) override
    {
        assert(idx < getNumVertices());
        return _vertices[idx];
    }

    // implements interface method
    const VertexInterface* getVertex(int idx) const override
    {
        assert(idx < getNumVertices());
        return _vertices[idx];
    }

 protected:
    VertexContainer _vertices;  //!< Vertex container

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template <class... Vertices>
class MixedEdge : public BaseMixedEdge
{
 public:
    using Ptr      = std::shared_ptr<MixedEdge>;
    using ConstPtr = std::shared_ptr<const MixedEdge>;
    using UPtr     = std::unique_ptr<MixedEdge>;

    //! Return number of vertices at compile-time
    static constexpr const int numVerticesCompileTime = sizeof...(Vertices);

    //! Typedef to represent the vertex container.
    using VertexContainer = std::array<VertexInterface*, numVerticesCompileTime>;

    // delete default constructor
    MixedEdge() = delete;

    /**
     * @brief Construct edge by providing connected vertices
     *
     * The order must match the order of template arguments provided to the class.
     * @warning Vertices must remain valid as long as the Edge is in use!
     */
    template <class... VerticesT>
    explicit MixedEdge(VerticesT&... args) : BaseMixedEdge(), _vertices({&args...})
    {
        // Check if types and number of template parameters match class template paremters for vertices
        static_assert(util::variadic_temp_equal<std::tuple<Vertices...>, std::tuple<VerticesT...>>::value,
                      "MixedEdge(): Number and types of vertices passed via the constructor does not match number of class template parameters");
    }  //!< Construct edge by passing all vertices references.

    // implements interface method
    int getObjectiveDimension() const override = 0;
    // implements interface method
    int getEqualityDimension() const override = 0;
    // implements interface method
    int getInequalityDimension() const override = 0;
    // implements interface method
    int getNumVertices() const override { return numVerticesCompileTime; }

    // implements interface method
    int verticesDimension() const override
    {
        int n = 0;
        for (VertexInterface* vertex : _vertices) n += vertex->getDimension();
        return n;
    }

    // implements interface method
    // bool isLinear() const override = 0;
    // implements interface method
    // bool providesJacobian() const override { return false; }

    // implements interface method
    void precompute() override                                                     = 0;
    void computeObjectiveValues(Eigen::Ref<Eigen::VectorXd> obj_values) override   = 0;
    void computeEqualityValues(Eigen::Ref<Eigen::VectorXd> eq_values) override     = 0;
    void computeInequalityValues(Eigen::Ref<Eigen::VectorXd> ineq_values) override = 0;

    // implements interface method
    VertexInterface* getVertexRaw(int idx) override
    {
        assert(idx < getNumVertices());
        return _vertices[idx];
    }

    // implements interface method
    const VertexInterface* getVertex(int idx) const override
    {
        assert(idx < getNumVertices());
        return _vertices[idx];
    }

 protected:
    const VertexContainer _vertices;  //!< Vertex container

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template <>
class MixedEdge<> : public BaseMixedEdge
{
 public:
    using Ptr      = std::shared_ptr<MixedEdge>;
    using ConstPtr = std::shared_ptr<const MixedEdge>;
    using UPtr     = std::unique_ptr<MixedEdge>;

    //! Typedef to represent the vertex container.
    using VertexContainer = std::vector<VertexInterface*>;

    // delete default constructor
    MixedEdge() : BaseMixedEdge() {}

    explicit MixedEdge(int num_vertices) : BaseMixedEdge() { resizeVertexContainer(num_vertices); }

    //! Set number \c n of vertices attached to this edge.
    void resizeVertexContainer(int n) { _vertices.resize(n); }

    void setVertex(int idx, VertexInterface& vertex)
    {
        assert(idx < (int)_vertices.size());
        _vertices[idx] = &vertex;
    }

    // implements interface method
    int getObjectiveDimension() const override = 0;
    // implements interface method
    int getEqualityDimension() const override = 0;
    // implements interface method
    int getInequalityDimension() const override = 0;
    // implements interface method
    int getNumVertices() const override { return (int)_vertices.size(); }

    // implements interface method
    int verticesDimension() const override
    {
        int n = 0;
        for (VertexInterface* vertex : _vertices) n += vertex->getDimension();
        return n;
    }

    // implements interface method
    // bool isLinear() const override = 0;
    // implements interface method
    // bool providesJacobian() const override { return false; }

    // implements interface method
    void precompute() override                                                     = 0;
    void computeObjectiveValues(Eigen::Ref<Eigen::VectorXd> obj_values) override   = 0;
    void computeEqualityValues(Eigen::Ref<Eigen::VectorXd> eq_values) override     = 0;
    void computeInequalityValues(Eigen::Ref<Eigen::VectorXd> ineq_values) override = 0;

    // implements interface method
    VertexInterface* getVertexRaw(int idx) override
    {
        assert(idx < getNumVertices());
        return _vertices[idx];
    }

    // implements interface method
    const VertexInterface* getVertex(int idx) const override
    {
        assert(idx < getNumVertices());
        return _vertices[idx];
    }

 protected:
    VertexContainer _vertices;  //!< Vertex container

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

}  // namespace corbo

#endif  // SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_EDGE_H_
