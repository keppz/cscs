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
                Gmsh.Option.SetNumber("General.Terminal", 1);
                Gmsh.Model.Add("t1");
                var lc = 1E-2;
                geo.AddPoint(0, 0, 0, lc, 1);
                geo.AddPoint(0.1, 0, 0, lc, 2);
                geo.AddPoint(0.1, 0.3, 0, lc, 3);
                var p4 = geo.AddPoint(0, 0.3, 0, lc);

                geo.AddLine(1, 2, 1);
                geo.AddLine(3, 2, 2);
                geo.AddLine(3, p4, 3);
                geo.AddLine(4, 1, p4);

                geo.AddCurveLoop(new int[] { 4, 1, -2, 3 }, 1);
                geo.AddPlaneSurface(new int[] { 1 }, 1);
                Gmsh.Model.AddPhysicalGroup(1, new int[] { 1, 2, 4 }, 5);
                var ps = Gmsh.Model.AddPhysicalGroup(2, new int[] { 1 });
                Gmsh.Model.SetPhysicalName(2, ps, "My surface");
                geo.Synchronize();
                Gmsh.Model.Mesh.Generate(2);

                Gmsh.Write("t1.msh");
                Gmsh.Fltk.Run();
                Gmsh.Finalize();
            }
            catch (Exception ex)
            {
                Status.Text = $"Error: {ex.Message}";
            }
        }
    }
}
