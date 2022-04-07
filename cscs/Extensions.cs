using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CustomExtensions
{
    public static class Extensions
    {
        public static IEnumerable<T> SliceCol<T>(this T[,] array, int col)
        {
            for (var i = 0; i < array.GetLength(0); i++)
            {
                yield return array[i, col];
            }
        }
        public static IEnumerable<T> SliceRow<T>(this T[,] array, int row)
        {
            for (var i = 0; i < array.GetLength(0); i++)
            {
                yield return array[row, i];
            }
        }
    }
}
