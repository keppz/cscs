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
using CustomExtensions;

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

                double lc = 1e-1;
                int nLoops;

                // Create a dummy array which is later passed to the function
                // call.
                List<double[,]> loops = new List<double[,]>();
                // Point coordinates in x, y, z
                double[,] loop1 = new double[4, 3]
                {
                    { 0.0d, 0.0d, 0.0d },
                    { 1.0d, 0.0d, 0.0d },
                    { 1.0d, 1.0d, 0.0d },
                    { 0.0d, 1.0d, 0.0d }
                };
                loops.Add(loop1);
                nLoops = loops.Count;
                int stag = 1;
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
                    int currentPointTag;
                    // Add lines
                    for (int iline = 0; iline < npts - 1; iline++)
                    {
                        currentPointTag = nlinesOffset + 1 + iline;
                        Gmsh.Model.Occ.AddLine(
                            currentPointTag,
                            currentPointTag + 1,
                            nlinesOffset + 1 + iline);
                        lineTags.Add(nlinesOffset + 1 + iline);
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
                    stag = Gmsh.Model.Occ.AddPlaneSurface(lltag, -1);
                    surfTags.Add(stag);
                    nptsOffset += npts;
                    nlinesOffset += npts;
                }

                // Write to model
                Gmsh.Model.Occ.Synchronize();

                // Create surfaces
                int surfTag = 1;
                (int, int)[] objectiveSurface = new (int, int)[1];
                (int, int)[] toolSurface = new (int, int)[1];
                (int, int)[] outDimTags;
                (int, int)[][] outDimTagsMap;
                if (nLoops > 1)
                {
                    for (var i = 1; i < nLoops; i++)
                    {
                        surfTag = nLoops + i;
                        objectiveSurface[0] = (2, surfTags[0]);
                        toolSurface[0] = (2, surfTags[i]);
                        Gmsh.Model.Occ.Cut(
                            objectiveSurface,
                            toolSurface,
                            out outDimTags,
                            out outDimTagsMap,
                            surfTag,
                            false, false);
                        surfTags[0] = surfTag;
                    }
                }
                else
                {
                    // only one surface, no cut performed
                    surfTag = stag;
                }

                // Syncronise: To build geometry
                Gmsh.Model.Occ.Synchronize();

                // remove unessary surfaces
                for (int i = 1; i < surfTags.Count; i++)
                {
                    (int, int)[] removeTags = { (2, i) };
                    Gmsh.Model.Occ.Remove(removeTags, false);
                    Gmsh.Model.Occ.Synchronize();
                }

                int[] pgtSurf = new int[nLoops];
                for (var i = 0; i < nLoops; i++)
                {
                    pgtSurf[i] = Gmsh.Model.AddPhysicalGroup(
                        1, lineLoopTags[i].ToArray(), -1);
                    Gmsh.Model.SetPhysicalName(
                        1, pgtSurf[i], $"LINE{i}");
                }

                int pgtEL;
                int[] physicalGroupTags = { surfTag };
                pgtEL = Gmsh.Model.AddPhysicalGroup(2, physicalGroupTags, -1);
                Gmsh.Model.SetPhysicalName(2, pgtEL, "ELvol");

                // Mesh the surface
                //
                // set mesh options

                // Delaunay
                Gmsh.Option.SetNumber("Mesh.Algorithm", 5);
                Gmsh.Option.SetNumber("Mesh.CharacteristicLengthMax", lc);
                Gmsh.Option.SetNumber("Mesh.ElementOrder", 2);

                // 1... msh 39...inp
                Gmsh.Option.SetNumber("Mesh.Format", 39);

                // number of smoothing steps
                Gmsh.Option.SetNumber("Mesh.Smoothing", 3);
                // save all node groups
                Gmsh.Option.SetNumber("Mesh.SaveGroupsOfNodes", 1);
                Gmsh.Option.SetNumber("Mesh.RecombineAll", 0);
                Gmsh.Option.SetNumber(
                    "Mesh.CharacteristicLengthExtendFromBoundary", 1);
                Gmsh.Option.SetNumber("Mesh.CharacteristicLengthMax", 1e22);
                Gmsh.Option.SetNumber("Mesh.CharacteristicLengthMin", lc);
                Gmsh.Option.SetNumber("Mesh.CharacteristicLengthFactor", 1);
                Gmsh.Option.SetNumber(
                    "Mesh.CharacteristicLengthFromPoints", 1);

                // Generate mesh
                Gmsh.Model.Mesh.Generate(1);
                Gmsh.Model.Mesh.Generate(2);

                long[] nodeTags;

                // stored in format [n1x, n1y, n1z, n2x, ...]
                double[] coordImport;
                double[] parametricCoord;

                // Get the nodes
                Gmsh.Model.Mesh.GetNodes(out nodeTags, out coordImport,
                    out parametricCoord, 2, -1, true, false);

                // reshape nodes
                int nNodes = nodeTags.Count();

                double[,] coords = new double[nNodes, 3];
                var perm = Enumerable.Range(0, nNodes).ToArray();
                // Sort array
                Array.Sort(nodeTags, perm);
                for (int row = 0; row < nNodes; row++)
                {
                    for (int col = 0; col < 3; col++)
                    {
                        coords[row, col] = coordImport[(perm[row] * 3) + col];
                    }
                }


                // Check if nodes are sorted
                for (int inode = 0; inode < nNodes - 1; inode++)
                {
                    if (nodeTags[inode] != nodeTags[inode + 1] - 1)
                    {
                        throw new InvalidOperationException(
                            "It is assumed that the nodes are ordered and " +
                            "contiguouse. Consult the Author for fixing this."
                            );
                    }
                }

                // Assume that nodes tags are sorted
                long minNodeTag = nodeTags.Min();

                // For debugging plot nodes
                var plt = new ScottPlot.Plot(800, 600);
                plt.AddScatter(
                    coords.SliceCol(0).ToArray(),
                    coords.SliceCol(1).ToArray(),
                    lineWidth: 0);
                plt.SaveFig("quickstart.png");

                // Get the Elements (There is only one type: 6 node tri)
                int[] elementTypes;

                // The element tags for each type [i][..]
                // The nodes for each type[i][..]. Stored in [e1n1, e1n2, ...]
                long[][] elementTags;

                long[][] nodeTagsInput2;
                Gmsh.Model.Mesh.GetElements(
                    out elementTypes, out elementTags,
                    out nodeTagsInput2, 2, -1);

                // Reshape elements array
                int nElm = elementTags[0].Count();
                const int nElNodes = 6;  // quadratic triangluars
                long[,] elements = new long[nElm, nElNodes];
                for (int ielm = 0; ielm < nElm; ielm++)
                {
                    for (int inode = 0; inode < nElNodes; inode++)
                    {
                        // Subtract the smallest nodeTag
                        elements[ielm, inode] =
                            nodeTagsInput2[0][(ielm * nElNodes) + inode] -
                            minNodeTag;
                    }
                }

                // Plot the elements
                int[] pntloop = { 0, 1, 2, 0 };
                for (int ielm = 0; ielm < nElm; ielm++)
                {
                    for (int iline = 0; iline < 3; iline++)
                    {
                        // Scottplot doesn't save data internally
                        // Therefore this data might get lost if it is
                        // out of scope
                        double[] lineX = new double[2];
                        double[] lineY = new double[2];
                        lineX[0] = coords[elements[ielm, pntloop[iline]], 0];
                        lineY[0] = coords[elements[ielm, pntloop[iline]], 1];
                        lineX[1] = coords[
                            elements[ielm, pntloop[iline + 1]], 0];
                        lineY[1] = coords[
                            elements[ielm, pntloop[iline + 1]], 1];
                        plt.AddScatter(lineX, lineY, lineWidth: 2);
                    }
                }

                double[] coordsTemp;
                long[] tagsSurf;

                plt.SaveFig("elements.png");
                // Get surfaces
                for (int isurf = 0; isurf < nLoops; isurf++)
                {
                    Gmsh.Model.Mesh.GetNodesForPhysicalGroup(
                        1, pgtSurf[isurf], out tagsSurf, out coordsTemp);
                }

                // Quit gmsh
                // Gmsh.Fltk.Run();
                Gmsh.Finalize();
            }
            catch (Exception ex)
            {
                Status.Text = $"Error: {ex.Message}";
            }
}
    }
}
