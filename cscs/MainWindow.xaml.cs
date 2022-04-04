using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using GmshNet;
using geo = GmshNet.Gmsh.Model.Geo;

namespace Cs
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml.
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
        }

        private void Button_Click(object sender, RoutedEventArgs e)
        {
            try
            {
                Gmsh.Initialize();
                Gmsh.Option.SetNumber("General.Terminal", 0);
                Gmsh.Model.Add("mesh");

                // Set surface visibility
                Gmsh.Option.SetNumber("Geometry.Surfaces", 1);
                Gmsh.Option.SetNumber("Geometry.PointNumbers", 1);
                Gmsh.Option.SetNumber(
                    "Geometry.OCCBooleanPreserveNumbering", 1);

                // Generate model
                int nptsOffset = 0;
                int nlinesOffset = 0;
                List<List<int>> lineLoopTags = new List<List<int>>();
                List<int> loopTags = new List<int>();
                List<int> urfTags = new List<int>();
                List<int> pointTags = new List<int>();
                List<int> lineTags = new List<int>();
                List<int> surfTags = new List<int>();

                double lc = 1e-3;
                int nLoops;

                // Create a dummy array which is later passed to the function
                // call.
                List<double[,]> loops = new List<double[,]>();
                double[,] loop1 = new double[4, 2]
                {
                    { 0, 0 },
                    { 1, 0 },
                    { 1, 1 },
                    { 1, 0 }
                };
                loops.Add(loop1);
                nLoops = loops.Count;
                for (var iloop = 0; iloop < nLoops; iloop++)
                {
                    // load data
                    var pts = loops[iloop];
                    var npts = pts.GetLength(0);

                    // Add nodes
                    for (int ipnt = 0; ipnt < npts; ipnt++)
                    {
                        pointTags.Add(nptsOffset + 1 + ipnt);
                        Gmsh.Model.Occ.AddPoint(
                            pts[ipnt, 0], pts[ipnt, 1], pts[ipnt, 2], lc,
                            pointTags[ipnt]);
                    }

                    // Add lines
                    for (int iline = 0; iline < npts - 1; iline++)
                    {
                        lineTags.Add(nlinesOffset + 1 + iline);
                        Gmsh.Model.Occ.AddLine(
                            pointTags[lineLoopTags[iline]],
                            pointTags[lineLoopTags[iline] + 1],
                            lineLoopTags[iline]);
                    }

                    // Add last line to close loop
                    lineTags.Add(npts + nlinesOffset);
                    Gmsh.Model.Occ.AddLine(pointTags[^1], pointTags[0],
                                           lineTags[^1]);

                    lineLoopTags.Add(lineTags);
                    int[] lltag = new int[1];
                    lltag[0] = Gmsh.Model.Occ.AddCurveLoop(
                        lineTags.ToArray(), -1);
                    loopTags.Add(lltag[0]);
                    var stag = Gmsh.Model.Occ.AddPlaneSurface(lltag, -1);
                    surfTags.Add(stag);
                    nptsOffset += npts;
                    nlinesOffset += npts;
                }

                // Write to model
                Gmsh.Model.Occ.Synchronize();
                if (nLoops > 1)
                {
                    for (var i = 1; i < nLoops; i++)
                    {
                        surf_tag = n_loops + i;
                        Gmsh.Model.Occ.Cut(
                            [(2, surf_tags[0])],
                            [(2, surf_tags[i])],
                            surf_tag, removeObject = False, removeTool = False)
                            surf_tags[0] = surf_tag
                    }
                }
                else
                {
                    // only one surface, no cut performed
                    surf_tag = stag
                }
                // Syncronise: To build geometry
                Gmsh.Model.Occ.Synchronize();
                // remove unessary surfaces
                for (i in surf_tag)
                {
                    Gmsh.Model.Occ.Remove([(2, i)], recursive = False);
                    Gmsh.Model.Occ.Synchronize();
                }

                for (var i = 0; i < n_loops; i++)
                {
                    pgt_ll = Gmsh.Model.AddPhysicalGroup(
                                        1, line_loop_tags[i], -1);
                    Gmsh.Model.SetPhysicalName(
                        1, pgt_ll, "NSETsurf%i" % (i));
                    pgt_ll = Gmsh.Model.AddPhysicalGroup(
                        1, line_loop_tags[i], -1);
                    Gmsh.Model.SetPhysicalName(
                        1, pgt_ll, "LINE%i" % i);
                }
                pgt_el = Gmsh.Model.AddPhysicalGroup(2, [surf_tag], -1);
                Gmsh.Model.SetPhysicalName(2, pgt_el, "ELvol");

                // Mesh the surface
                //
                // set mesh options

                // Delaunay
                Gmsh.Option.SetNumber("Mesh.Algorithm", 5);
                Gmsh.Option.SetNumber("Mesh.CharacteristicLengthMax", lc);
                Gmsh.Option.SetNumber("Mesh.ElementOrder", 2)
                Gmsh.Option.SetNumber("Mesh.Format", 39)  # 1... msh 39...inp
                Gmsh.Option.SetNumber("Mesh.Smoothing", 3)  # number of smoothing steps
                Gmsh.Option.SetNumber("Mesh.SaveGroupsOfNodes", 1)  # save all node groups
                Gmsh.Option.SetNumber("Mesh.RecombineAll", 0)
                Gmsh.Option.SetNumber("Mesh.CharacteristicLengthExtendFromBoundary", 1)
                Gmsh.Option.SetNumber("Mesh.CharacteristicLengthMax", 1e22)
                Gmsh.Option.SetNumber("Mesh.CharacteristicLengthMin", lc)
                Gmsh.Option.SetNumber("Mesh.CharacteristicLengthFactor", 1)
                Gmsh.Option.SetNumber("Mesh.CharacteristicLengthFromPoints", 1)

                // Generate mesh
                Gmsh.model.mesh.generate(1);
                Gmsh.model.mesh.generate(2);
                // write
                // Quit gmsh
                Gmsh.Fltk.Run();
                Gmsh.Finalize();

                int dim = 2;
                int tag = -1;
                bool includeBoundary = false;
                bool returnParametricCorrd = false;
                // Utilize the routines getNodes(dim, tag, includeBoundary, returnParametricCorrd),
                // std::vector<std::size_t> nodeTags;
                // std::vector<double> nodeCoords, nodeParams;
                // gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, dim, tag);

                // Get the mesh elements for the entity (dim, tag):
                // std::vector<int> elemTypes;
                // std::vector<std::vector<std::size_t>> elemTags, elemNodeTags;
                // gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, dim, tag);

                Gmsh.Model.Mesh.GetNodes(dim, tag, includeBoundary, returnParametricCoord);
            }
            catch (Exception ex)
            {
                Status.Text = $"Error: {ex.Message}";
            }
        }
    }
}
