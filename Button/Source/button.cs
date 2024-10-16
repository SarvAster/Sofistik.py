using System.Diagnostics;
using System.Windows.Forms;

public class RibbonAddin
{
    public void OnThermalButtonClick(Office.IRibbonControl control)
    {
        // Path to the Python interpreter
        string pythonExe = @"C:\Path\To\python.exe";

        // Path to the Python script
        string pythonScript = @"C:\Path\To\Scripts\modify_dat_file.py";

        // Path to the .dat file
        string datFilePath = @"C:\Path\To\your_project.dat";

        ProcessStartInfo start = new ProcessStartInfo();
        start.FileName = pythonExe;
        start.Arguments = string.Format("\"{0}\" \"{1}\"", pythonScript, datFilePath);
        start.UseShellExecute = false;
        start.RedirectStandardOutput = true;
        start.RedirectStandardError = true;
        start.CreateNoWindow = true;

        using (Process process = Process.Start(start))
        {
            using (System.IO.StreamReader reader = process.StandardOutput)
            {
                string result = reader.ReadToEnd();
                MessageBox.Show(result);  // Show any output or errors
            }
        }
    }
}
