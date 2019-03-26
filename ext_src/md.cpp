#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/split.hpp>


namespace py = pybind11;
using namespace std;
using namespace boost;
typedef vector< string > split_vector_type;
using boost::lexical_cast;


std::vector<int> merge_main_output(const std::vector<std::string> &folders, 
     const std::string& merged_filename, bool verbose)
{
    std::string filename = "dftb.out";
    vector<int> last_step_nos;
    if (verbose)
    {
        std::cout << "Combining */" << filename << " to " << merged_filename << std::endl;
    }
    
    ofstream fout( merged_filename );
    bool write_header = false;
    for (size_t i=0; i<folders.size(); ++i)
    {
        string filein = (folders[i] + "/" + filename);
        if (verbose)
            cout <<  "  Loading " << filein << endl;
        ifstream fin( filein.c_str() );
        stringstream buffer; 
        string str;

        int max_step = 0;
        int min_step = 0;
        bool first = true;

        while (getline(fin, str))
        {
            if (starts_with(str, "  *** Start molecular dynamics ***"))
            {
                if (write_header == false)
                {
                    write_header = true;
                    buffer << str << endl;
                    fout << buffer.str() << endl;
                }
                buffer.str("");
            }
            else if (starts_with(str, " *** AT T="))
            {
                vector<string> t_lines;
                t_lines.push_back(str);             

                split_vector_type SplitVec; 
                split( SplitVec, str, is_any_of(" "), token_compress_on );
                int step_no = lexical_cast<int>(SplitVec[10]) ;
                if (first)
                {
                    min_step = step_no;
                    first = false;
                }
                try
                {
                    for (int i=0; i<4;++i)
                    {
                        getline(fin, str);
                        t_lines.push_back(str);
                    }
                    if (step_no > max_step)
                    {
                        max_step = step_no;
                    }                        
                    for (const auto & s : t_lines)
                    {
                        buffer << s << endl;
                    }
                    bool to_write = true;
                    if (!last_step_nos.empty())
                    {
                        if (step_no <= last_step_nos.back())
                        {
                            to_write = false;
                        }
                    }
                    if (to_write)
                    {
                        fout << buffer.str() << endl;
                        buffer.str("");
                    }
                    else
                    {
                        buffer.str("");
                    }
                }
                catch (std::ifstream::failure e) {
                    break;
                }
            }
            else
            {
                buffer << str << endl;
            }
        }
        fout << buffer.str() << endl;
        if (verbose)
            cout << min_step << max_step << endl;
        last_step_nos.push_back(max_step);
    }
    return last_step_nos;
}

void merge_datafile(const std::vector<std::string> &folders, const vector<int>& last_step_nos,
    const std::string& filename, const std::string& merged_filename, bool verbose)
{
    string str;
    if (verbose)
        cout << "Combining */" << filename << " to merged_" << merged_filename  << endl;

    ofstream fout( ("merged_"+merged_filename).c_str() );
    for (size_t i=0; i<folders.size(); ++i)
    {
        string filein = (folders[i] + "/" + filename );
        if (verbose)
            cout <<  "  Loading " << filein << endl;

        ifstream fin( filein.c_str() );

        bool first_line = true;

        int nat = 0;
        while (getline(fin, str))
        {
            string natline = str;
            // fout << str << endl;  //nat
            
            if (first_line)
            {
                first_line = false;
                split_vector_type SplitVec; 
                trim(str);
                split( SplitVec, str, is_any_of(" "), token_compress_on );    
                nat = lexical_cast<int>(SplitVec[0]) ;                
            }          
                        
            getline(fin, str);
            string title_line = str;
            // fout << str << endl;  //title

            bool to_write = true;
            
            split_vector_type SplitVec; 
            split( SplitVec, str, is_any_of(" "), token_compress_on );
            int step_no = lexical_cast<int>(SplitVec[10]) ;

            // cout << step_no << endl;
            if (i==0)
            {
                if (step_no > last_step_nos[i])
                {
                    break;
                }
            }
            else if ( step_no <= last_step_nos[i-1] )
            {
                to_write = false;
            }
            if (to_write)
            {
                fout << natline << endl;
                fout << title_line << endl;
            }
            for (int j=0; j<nat;++j)
            {
                getline(fin, str);
                if (to_write)
                {
                    fout << str << endl;
                }
            }   
        }
    }
}

void stride_xyz(const std::string& input, const std::string& output, int size)
{
    ifstream fin(input);
    ofstream fout(output);

    int stride_size = size;
    string str;
    int ind = 0;
    while (getline(fin, str))
    {
        if (ind % stride_size == 0)
        {
            trim(str);
            fout << str << endl;
            int nat = boost::lexical_cast<int>(str);
            getline(fin, str);
            fout << str << endl;
            for (int i=0; i<nat; ++i)
            {
                getline(fin, str);
                fout << str << endl;
            }
        }
        else
        {
            trim(str);
            int nat = boost::lexical_cast<int>(str);
            getline(fin, str);
            for (int i=0; i<nat; ++i)
            {
                getline(fin, str);
            }
        }
        ind += 1;
    }
}

PYBIND11_MODULE(_dcdftbmd_tools_cpp_extension_md, m) {
    m.doc() = R"pbdoc(
        
    )pbdoc";

    m.def("merge_main_output", &merge_main_output, R"pbdoc(
        Merge the main output (dftb.out) of DCDFTBMD code.

    )pbdoc");

    m.def("merge_datafile", &merge_datafile, R"pbdoc(
        Merge a specific data file (e.g., trajectory file) of DCDFTBMD.

        
    )pbdoc");

    m.def("stride_xyz", &stride_xyz, R"pbdoc(
        Stride trajectory file of DCDFTBMD.

        
    )pbdoc");

    

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
