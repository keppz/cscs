using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Cs
{
    internal class FemSolver
    {

        //
        // Meshing function
        //
        public void ReadInput(int boundaryCurves)
        {
            throw new NotSupportedException();
        }

        // ====================================================================
        // Transformation function
        // ====================================================================
        private double trans(double t, double a, double b)
        {
            // Transform the line element border.
            return (a * t) + b;
        }                                                                      

//def solve_cross_section() :
//    """Mainprogram routine."""
//    # Some predefined constants
//    filename = './mesh/out.inp'  # input file
//    control = {'plot': False
//    }
//# Read the input file
//    nodes, elements, surfaces = read_input(filename)
//    nel_nodes = len(elements[0, 1:])
//    n_nodes = len(nodes)
//    # inner boundaries
//    n_loops = len(surfaces)
//    # get element order
//    if nel_nodes == 3:
//        el_order = 1
//        n_edge = 2
//        n_int_1d = 1
//        n_int = 1
//    elif nel_nodes == 6:
//        el_order = 2
//        n_edge = 3
//        n_int_1d = 2
//        n_int = 3

//    # -------------------------------------------------------------------------
//    # Get integration points
//    gaup, weights = g_2d(el_order)
//    gaup_1d, weights_1d = g_1d(2)
//    gaup_area, weights_area = g_2d(3)
//    n_int_area = 7
//    n_int_1d = 2
//    # initalise array
//    c_nodes = []
//    for mm in range(n_loops) :
//        c_nodes.append(np.zeros((n_edge* len(surfaces[mm]), 1), dtype=int))

//    # =========================================================================
//    # FEM analysis
//    # =========================================================================
//    # -------------------------------------------------------------------------
//    # Allocate global stiffnes matrix
//    # -------------------------------------------------------------------------
//    r, c = elements.shape
//    n_elm = r
//    r, c = nodes.shape
//    n_node = r

//    f_global = np.zeros((n_node, 1))
//    # indexing arrays
//    row = np.zeros((n_elm* nel_nodes * nel_nodes, 1))
//    col = np.zeros((n_elm* nel_nodes * nel_nodes, 1))
//    # stiffness matrice array
//    val = np.zeros((n_elm* nel_nodes * nel_nodes, 1))
//    # -------------------------------------------------------------------------
//    # Loop over elements, compute the local stiffness matrix and assemble them
//    # in to the global stiffness matrix
//    # -------------------------------------------------------------------------
//    u = np.zeros((n_nodes, n_loops))
//    for i in range(n_elm) :
//        # ---------------------------------------------------------------------
//        # Get nodal coordinates of the element
//        # remove z coordinate: Your mesh must lie in the x-y plane!
//        el_nodes = elements[i, 1:]
//        el_node_coord = nodes[el_nodes - 1, 1:3]
//        k_el = np.zeros((nel_nodes, nel_nodes))
//        f_el = np.zeros((nel_nodes, 1))
//        # numerical integration
//        for g in range(n_int) :
//            gp = gaup[g, 1:]  # current gauss point
//            Nsp = dNi(gp, el_order)

//            JKinv = [[np.dot(el_node_coord[:, 0], Nsp[:, 0]),
//                      np.dot(el_node_coord[:, 0], Nsp[:, 1])],
//                     [np.dot(el_node_coord[:, 1], Nsp[:, 0]),
//                      np.dot(el_node_coord[:, 1], Nsp[:, 1])]]
//            JK = np.linalg.inv(JKinv)
//            detJKinv = np.linalg.det(JKinv)
//            # -----------------------------------------------------------------
//            # element vector
//            # -----------------------------------------------------------------
//            bb = np.matmul(Nsp, JK)
//            k_el = k_el + weights[g, 1] * detJKinv * np.matmul(bb, bb.T)
//            # -----------------------------------------------------------------
//            # element vector
//            # -----------------------------------------------------------------
//            NKs = Ni(gp, el_order)
//            f_el = f_el + weights[g, 1] * detJKinv * (-1) * NKs

//        # ---------------------------------------------------------------------
//        # Detect if element is part of surface
//        # ---------------------------------------------------------------------
//        for m in range(len(surfaces)) :
//            ind = np.where(surfaces[m][:, 0] == i + 1)
//            ind = ind[0]
//            if ind.size == 0:
//                continue

//            sur = surfaces[m][ind, 1]
//            for nn in range(len(sur)) :
//                # -------------------------------------------------------------
//                # get nodes that lie on the surface
//                if sur[nn] == 1:
//                    edg_nodes = [0, 1, 3]
//    elif sur[nn] == 2:
//                    edg_nodes = [1, 2, 4]
//    elif sur[nn] == 3:
//                    edg_nodes = [2, 0, 5]

//# -------------------------------------------------------------
//# write node node number into boundary array
//    c_nodes[m][n_edge * (ind - 1)] = el_nodes[edg_nodes[0]] - 1
//                c_nodes[m][n_edge * (ind - 1) + 1] = el_nodes[edg_nodes[1]] - 1
//                if el_order == 2:
//                    c_nodes[m][n_edge * (ind - 1) + 2] = \
//                                                     el_nodes[edg_nodes[2]] - 1

//        # ---------------------------------------------------------------------
//        # Assemble into global matrix
//        # ---------------------------------------------------------------------
//        index_global = el_nodes - 1
//        for l in range(nel_nodes) :
//            f_global[index_global[l]] = f_global[index_global[l]] + f_el[l]
//            for k in range(nel_nodes) :
//                row[i * nel_nodes * *2 + nel_nodes * l + k] = index_global[l]
//                col[i * nel_nodes * *2 + nel_nodes * l + k] = index_global[k]
//                val[i * nel_nodes * *2 + nel_nodes * l + k] = k_el[l, k]

//# =========================================================================
//# Include boundary values by modifing value array
//# =========================================================================
//    cond_factor = val.max()

//    for i in range(len(c_nodes)) :
//        for j in range(len(c_nodes[i])) :
//            node_nr = c_nodes[i][j]
//            val[row == node_nr] = 0
//            indexes = np.where((col == node_nr) & (row == node_nr))
//            # set only one entry since they are summed later anyways
//            val[indexes[0][0]] = cond_factor

//# =========================================================================
//# Write values into sparse global matrix
//# =========================================================================
//# double indexed entries are summed up
//    k_global = scipy.sparse.coo_matrix((val.flatten(), (row.flatten(),
//                                        col.flatten())),
//                                       shape = (n_node, n_node))
//    # =========================================================================
//    # Boundary
//    # =========================================================================
//    # Remove non unique nodes
//    for i in range(len(c_nodes)) :
//        c_nodes[i] = np.unique(c_nodes[i])

//    # -------------------------------------------------------------------------
//    # Modify rhs
//    for mm in range(n_loops) :
//        mag = np.zeros((n_loops, 1))
//        if mm > 0:
//            mag[mm] = 1

//        # ---------------------------------------------------------------------
//        # Boundary force vector
//        for i in range(n_loops) :
//            for j in range(len(c_nodes[i])) :
//                f_global[c_nodes[i][j]] = mag[i] * cond_factor

//        # =====================================================================
//# Solve equation system
//# =====================================================================
//    u[:, mm] = scipy.sparse.linalg.spsolve(k_global.tocsr(), f_global)
//    # =========================================================================
//    # Linear superposition principle
//    # =========================================================================
//    # Build a lhs matrich
//    # lhs{1}(:)*u(:,1) + lhs{1}(:)*u(:,2) = rhs(1) # first loop
//    # lhs{2}(:)*u(:,1) + lhs{2}(:)*u(:,2) = rhs(2) # second loop
//    # ...
//    lhs = []
//    psi_loop = []
//    for i in range(n_loops) :
//        nsurf = len(surfaces[i])
//        lhs.append(np.zeros((1, nel_nodes* nsurf)))
//        psi_loop.append(np.zeros((nel_nodes* nsurf, n_loops)))

//    rhs = np.zeros((n_loops, 1))
//    for i in range(n_elm) :
//        # ---------------------------------------------------------------------
//        # Get nodal coordinates of the element
//        # remove z coordinate: Your mesh must lie in the x-y plane!
//        el_nodes = elements[i, 1:]
//        el_node_coord = nodes[el_nodes - 1, 1:3]
//        # ---------------------------------------------------------------------
//        # Detect if edges are part of a surface
//        # only for the inner surfaces
//        for m in range(n_loops) :
//            ind = np.where(surfaces[m][:, 0] == i + 1)
//            ind = ind[0]
//            if ind.size == 0:
//                continue

//            sur = surfaces[m][ind, 1]
//            for nn in range(len(sur)) :
//                sind = ind[nn]
//                psi_loop[m][(sind * nel_nodes):
//                            (sind * nel_nodes + nel_nodes), :] = \
//                                                             u[el_nodes - 1, :]
//                # -------------------------------------------------------------
//                # Get the coordinates of the surface
//                if sur[nn] == 1:
//                    xi_1 = np.array([0, 0])
//                    xi_2 = np.array([1, 0])
//                elif sur[nn] == 2:
//                    xi_1 = np.array([1, 0])
//                    xi_2 = np.array([0, 1])
//                elif sur[nn] == 3:
//                    xi_1 = np.array([0, 1])
//                    xi_2 = np.array([0, 0])

//                # -------------------------------------------------------------
//                # Line integration
//                a = (xi_2 - xi_1)
//                b = xi_1
//                for g in range(n_int_1d) :
//                    # Get 2D positon of Gauss point
//                    gp = trans(gaup_1d[g, 1], a, b)
//                    # ---------------------------------------------------------
//                    # Geometric information at this position
//                    Ng = Ni(gp, el_order)
//                    # tangent
//                    dNt = dNti(gaup_1d[g, 1], a, b, el_order)
//                    #
//                    x = np.matmul(el_node_coord.T, Ng)
//                    dxdt = np.matmul(el_node_coord.T, dNt)
//                    normal = np.flipud(dxdt)
//                    rv = dxdt / np.linalg.norm(dxdt)
//                    normal[1] = -normal[1]
//                    normal = normal / np.linalg.norm(normal)
//                    Nsp = dNi(gp, el_order)
//                    JKinv = [[np.dot(el_node_coord[:, 0], Nsp[:, 0]),
//                              np.dot(el_node_coord[:, 0], Nsp[:, 1])],
//                             [np.dot(el_node_coord[:, 1], Nsp[:, 0]),
//                              np.dot(el_node_coord[:, 1], Nsp[:, 1])]]
//                    JK = np.linalg.inv(JKinv)
//# ---------------------------------------------------------
//# Compute the rhs and lhs
//# rhs gives 2 the area of the loop ( good to check things)
//                    rhs[m] = rhs[m] + weights_1d[g, 1] * (x[0] * rv[1]
//                                    - x[1] * rv[0]) * np.linalg.norm(dxdt)
//                    for k in range(nel_nodes) :
//                        lhs[m][0, sind * nel_nodes + k] = \
//                            lhs[m][0, sind * nel_nodes + k] \
//                            + weights_1d[g, 1] * 2 \
//                            * np.dot(np.matmul(Nsp[k, :], JK), normal) \
//                            * np.linalg.norm(dxdt)

//    # -------------------------------------------------------------------------
//# Compute superposition factors
//    lhs_mat = np.zeros((n_loops, n_loops))
//    for i in range(n_loops) :
//        for j in range(n_loops) :
//            lhs_mat[i, j] = np.dot(lhs[i][0, :], psi_loop[i][:, j])

//    mf = np.linalg.solve(lhs_mat, rhs)
//    # =========================================================================
//    # Final solution
//    # =========================================================================
//    psi = np.zeros((n_nodes, 1))
//    for i in range(n_loops) :
//        psi = psi + mf[i] * u[:, i, None]

//    # =========================================================================
//# Compute area integral for torsional inertia moment
//# =========================================================================
//    ij = np.zeros((n_node, 1))
//    for i in range(n_elm) :
//        # ---------------------------------------------------------------------
//        # Get nodal coordinates of the element
//        # remove z coordinate: Your mesh must lie in the x-y plane!
//        el_nodes = elements[i, 1:]
//        el_node_coord = nodes[el_nodes - 1, 1:3]
//        IKs = np.zeros((nel_nodes, 1))
//        TKinv = np.zeros((2, 1))
//        for g in range(n_int) :
//            gp = gaup[g, 1:]  # current gauss point
//            NKs = Ni(gp, el_order)
//            Nsp = dNi(gp, el_order)
//            JKinv = [[np.dot(el_node_coord[:, 0], Nsp[:, 0]),
//                      np.dot(el_node_coord[:, 0], Nsp[:, 1])],
//                     [np.dot(el_node_coord[:, 1], Nsp[:, 0]),
//                      np.dot(el_node_coord[:, 1], Nsp[:, 1])]]
//            detJKinv = np.linalg.det(JKinv)
//            if n_loops == 1:
//                IKs = IKs + (-4) * weights[g, 1] * detJKinv * NKs
//            else:
//                JK = np.linalg.inv(JKinv)
//                TKinv[1, 0] = np.dot(el_node_coord[:, 1], NKs)
//                TKinv[0, 0] = np.dot(el_node_coord[:, 0], NKs)
//                IKs = IKs + 2 * weights[g, 1] * np.matmul(np.matmul(Nsp, JK),
//                                                          TKinv) * detJKinv

//        # assemble
//    index_global = el_nodes - 1
//        for l in range(nel_nodes) :
//            ij[index_global[l]] = ij[index_global[l]] + IKs[l]

//# -------------------------------------------------------------------------
//# Compute area moments
//    ixx = 0.
//    iyy = 0.
//    ixy = 0.
//    qx = 0.
//    qy = 0.
//    area = 0.
//    for i in range(n_elm) :
//        el_nodes = elements[i, 1:]
//        el_node_coord = nodes[el_nodes - 1, 1:3]
//        IKs = np.zeros((nel_nodes, 1))
//        TKinv = np.zeros((2, 1))
//        for g in range(n_int_area) :
//            gp = gaup_area[g, 1:]  # current gauss point
//            NKs = Ni(gp, el_order)
//            Nsp = dNi(gp, el_order)
//            JKinv = [[np.dot(el_node_coord[:, 0], Nsp[:, 0]),
//                      np.dot(el_node_coord[:, 0], Nsp[:, 1])],
//                     [np.dot(el_node_coord[:, 1], Nsp[:, 0]),
//                      np.dot(el_node_coord[:, 1], Nsp[:, 1])]]
//            detJKinv = np.linalg.det(JKinv)

//            x_pos = np.dot(el_node_coord[:, 0], NKs)
//            y_pos = np.dot(el_node_coord[:, 1], NKs)
//            darea_el = weights_area[g, 1] * detJKinv
//            area += darea_el
//            ixx += darea_el* y_pos ** 2
//            iyy += darea_el* x_pos ** 2
//            ixy += darea_el* x_pos * y_pos
//            qx += darea_el* y_pos
//            qy += darea_el* x_pos

//    # -------------------------------------------------------------------------
//# compute torsional moment of inertia
//    it = np.matmul(psi.T, ij)
//# return to calc_torsional_constant and then to gui, so there plotted
//    cross_section_props = {
//            'it': it, 'ixx': ixx, 'iyy': iyy, 'ixy': ixy,
//            'area': area, 'qx': qx, 'qy': qy}
//    return cross_section_props
//# =========================================================================
//# end of program
//# =========================================================================

//# =============================================================================
//# Gauss points
//# =============================================================================
//def g_1d(order) :
//    int_points = np.array([[1, 0.5 * (1 - np.sqrt(1 / 3))],
//                           [2, 0.5 * (1 + np.sqrt(1 / 3))]])
//    weights = np.array([[1, 0.5],
//                        [2, 0.5]])
//    return int_points, weights


//def g_2d(order) :
//    if order == 1:
//        int_points = np.array([[1, 1 / 3, 1 / 3], ])  # [ gauss point number, x, y]
//        weights = np.array([[1, 1 / 2], ])  # [ gauss point number, weight]

//    if order == 2:
//        int_points = np.array([[1, 0.1666666666666667, 0.1666666666666667],
//                               [2, 0.6666666666666667, 0.1666666666666667],
//                               [3, 0.1666666666666667, 0.6666666666666667]])
//        weights = np.array([[1, 0.3333333333333333 / 2],  # factor 0.5
//                            [2, 0.3333333333333333 / 2],
//                            [3, 0.3333333333333333 / 2]])
//    if order == 3:
//        int_points = np.array([[1, 0.1012865073235, 0.1012865073235],
//                               [2, 0.7974269853531, 0.1012865073235],
//                               [3, 0.1012865073235, 0.7974269853531],
//                               [4, 0.4701420641051, 0.0597158717898],
//                               [5, 0.4701420641051, 0.4701420641051],
//                               [6, 0.0597158717898, 0.4701420641051],
//                               [7, 0.3333333333333, 0.3333333333333]])
//        weights = np.array([[1, 0.1259391805448 / 2],
//                            [2, 0.1259391805448 / 2],
//                            [3, 0.1259391805448 / 2],
//                            [4, 0.1323941527885 / 2],
//                            [5, 0.1323941527885 / 2],
//                            [6, 0.1323941527885 / 2],
//                            [7, 0.2250000000000 / 2]])

//    return int_points, weights

//# =============================================================================
//# Triangle Elements
//# =============================================================================
//# linear
//#
//# X2       n3 o                                       | xi_2
//# |          /\          * ... integration points     |\
//# |         /  \         o ... nodes                  | \
//# |        /    \                                     |  \
//# |       /      \  n2                                |   \
//# |   n1 o________o___                                | 1* \
//# |____________________>X1                            |_____\_____>xi_1
//#
//# ansatz functions
//# gaup = [xi_1, xi_2]^T
//#
//# quadratic
//#
//# X2       n3 o                                       | xi_2
//# |          /\          * ... integration points     |\
//# |         /  \         o ... nodes                  | \
//# |    n6  o    o n5                                  | *3
//# |       /      \  n2                                |   \
//# |   n1 o____o____o___                               |1* *2
//# |___________n4_________>X1                          |_____\_____> xi_1
//#
//# ansatz functions
//# gaup = [xi_1, xi_2]^T


//def Ni(g, el_order) :
//    if el_order == 1:
//        return np.array([[1 - g[0] - g[1]],
//                         [g[0]],
//                         [g[1]]])
//    elif el_order == 2:
//        return np.array(
//            [[1 - g[0] - g[1] - 0.5 * (4 * g[0] * (1 - g[0] - g[1]))        # 1
//              - 0.5 * (4 * g[1] * (1 - g[0] - g[1]))],
//             [g[0] - 0.5 * (4 * g[0] * (1 - g[0] - g[1]))                   # 2
//              - 0.5 * 4 * g[0] * g[1]],
//             [g[1] - 0.5 * (4 * g[1] * (1 - g[0] - g[1]))                   # 3
//              - 0.5 * 4 * g[0] * g[1]],
//             [4 * g[0] * (1 - g[0] - g[1])],                                # 4
//             [4 * g[0] * g[1]],                                             # 5
//             [4 * g[1] * (1 - g[0] - g[1])]])                               # 6


//# =============================================================================
//# Derivatives of anatz functions
//# =============================================================================
//def dNi(g, el_order):
//    if el_order == 1:

//        return np.array([[-1, -1],
//                         [1, 0],
//                         [0, 1]])
//        # integration points in unit element
//    elif el_order == 2:
//        # Derivatives of anatz functions
//        return np.array([
//            [2.0 * g[1] - 2.0 * (-g[1] - g[0] + 1) + 2.0 * g[0] - 1,
//             2.0 * g[1] - 2.0 * (-g[1] - g[0] + 1) + 2.0 * g[0] - 1],       # 1
//            [-2.0 * g[1] - 2.0 * (-g[1] - g[0] + 1) + 2.0 * g[0] + 1, 0],   # 2
//            [0, 2.0 * g[1] - 2.0 * (-g[1] - g[0] + 1) - 2.0 * g[0] + 1],    # 3
//            [4 * (-g[1] - g[0] + 1) - 4 * g[0], -4 * g[0]],                   # 4
//            [4 * g[1], 4 * g[0]],                                           # 5
//            [-4 * g[1], 4 * (-g[1] - g[0] + 1) - 4 * g[1]]])                # 6


//# =============================================================================
//# line element derivation
//# =============================================================================
//def dNti(t, A, B, el_order):
//    if el_order == 1:
//        return np.array([[-A[1] - A[0]],
//                         [+A[0]],
//                         [+A[1]]])
//    elif el_order == 2:
//        return np.array([
//            [-2.0 * (-A[1] - A[0]) * (A[1] * t + B[1])
//             - 2.0 * A[1] * (-A[1] * t - A[0] * t - B[1] - B[0] + 1)
//             - 2.0 * A[0] * (-A[1] * t - A[0] * t - B[1] - B[0] + 1)
//             - 2.0 * (-A[1] - A[0]) * (A[0] * t + B[0]) - A[1] - A[0]],     # 1
//            [-2.0 * A[0] * (A[1] * t + B[1]) \
//             -2.0 * A[0] * (-A[1] * t - A[0] * t - B[1] - B[0] + 1)
//             - 2.0 * A[1] * (A[0] * t + B[0]) \
//             -2.0 * (-A[1] - A[0]) * (A[0] * t + B[0]) + A[0]],            # 2
//            [-2.0 * (-A[1] - A[0]) * (A[1] * t + B[1])
//             - 2.0 * A[0] * (A[1] * t + B[1])
//             - 2.0 * A[1] * (-A[1] * t - A[0] * t - B[1] - B[0] + 1)
//             - 2.0 * A[1] * (A[0] * t + B[0]) + A[1]],                      # 3
//            [+4 * A[0] * (-A[1] * t - A[0] * t - B[1] - B[0] + 1)
//             + 4 * (-A[1] - A[0]) * (A[0] * t + B[0])],                     # 4
//            [+4 * A[0] * (A[1] * t + B[1])
//             + 4 * A[1] * (A[0] * t + B[0])],                               # 5
//            [+4 * (-A[1] - A[0]) * (A[1] * t + B[1])
//             + 4 * A[1] * (-A[1] * t - A[0] * t - B[1] - B[0] + 1)]])      # 6
    }
}
