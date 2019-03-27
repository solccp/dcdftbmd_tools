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
#include <boost/format.hpp>



namespace py = pybind11;
using namespace std;
using namespace boost;
typedef vector< string > split_vector_type;
using boost::lexical_cast;
using boost::format;
using boost::io::group;

// parse and merge the main dcdftbmd output
// args:
//   folders: folders to be merged, must be in order.
//   input_filename: the name of the dcdftbmd output filename, should be dftb.out
//   verbose: print additional information
//   merged_filename: merged dcdftbmd main output filename
//   write_merged: write merged_filename
//   merged_data_filename: filename for merged parsed data
//   write_parsed_data: write merged parsed data
std::vector<int> merge_main_output(const std::vector<std::string> &folders, 
     const std::string& input_filename, bool verbose, int stride_size,
     const std::string& merged_filename, bool write_merged, 
     const std::string& merged_data_filename, bool write_parsed_data){
    vector<int> last_step_nos;
    
    ofstream fout;
    if (write_merged)
    {
        fout.open(merged_filename);
    }

    ofstream fout_data;
    if (write_parsed_data)
    {
        fout_data.open(merged_data_filename);
        fout_data << format("%1$-10s %2$-7s %3$-8s %4$-14s %5$-20s %6$-20s %7$-20s") % "#Time(fs)" % "Step" %
                    "#Sub.Sys" % "Temperature(K)" % "Pot. Energy(H)" % "Kinetic Energy(H)" % "MD Energy(H)" << endl; 
    }

    bool write_header = false;
    for (size_t i=0; i<folders.size(); ++i)
    {
        string filein = (folders[i] + "/" + input_filename);
        if (verbose)
        {
            cerr <<  "  Loading " << filein << ": ";
        }
            
        ifstream fin( filein.c_str() );
        stringstream buffer; 
        string str;

        int max_step = 0;
        int min_step = 0;
        bool first = true;

        //data section
        int no_subsystems;
        int step_no;
        double step_time;
        double step_temperature;
        double step_kinEne;
        double step_mdEne;
        double step_potEne;

        while (getline(fin, str))
        {
            if (starts_with(str, "  *** Start molecular dynamics ***"))
            {
                if (write_header == false)
                {
                    write_header = true;
                    buffer << str << endl;
                    if (write_merged)
                        fout << buffer.str() << endl;
                }
                buffer.str("");
            }
            else if (starts_with(str,"  ***  Number of subsystems ="))
            {
                split_vector_type SplitVec; 
                split( SplitVec, str, is_any_of(" "), token_compress_on );
                no_subsystems = lexical_cast<int>(SplitVec[6]) ;
                buffer << str << endl;
            }
            else if (starts_with(str, "    Final") && contains(str, "iterations"))
            {
                split_vector_type SplitVec; 
                split( SplitVec, str, is_any_of(" "), token_compress_on );
                step_potEne = lexical_cast<double>(SplitVec[5]) ;
                buffer << str << endl;
            }
            else if (starts_with(str, " *** AT T="))
            {
                vector<string> t_lines;
                t_lines.push_back(str);             

                split_vector_type SplitVec; 
                split( SplitVec, str, is_any_of(" "), token_compress_on );
                step_no = lexical_cast<int>(SplitVec[10]) ;
                step_time = lexical_cast<double>(SplitVec[4]) ;

                if (first)
                {
                    min_step = step_no;
                    first = false;
                    if (last_step_nos.size() > 0)
                    {
                        if (step_no < last_step_nos.back())
                        {
                            throw runtime_error("The order of folder to merge is wrong.");
                        }
                    }
                }

                
                try
                {
                    for (int i=0; i<4;++i)
                    {
                        getline(fin, str);
                        split_vector_type SplitVec; 
                        split( SplitVec, str, is_any_of(" "), token_compress_on );
                        switch (i)
                        {
                            case 1:
                                step_temperature = lexical_cast<double>(SplitVec[3]) ;
                                break;
                            case 2:
                                step_kinEne = lexical_cast<double>(SplitVec[4]) ;
                                break;
                            case 3:
                                step_mdEne = lexical_cast<double>(SplitVec[5]) ;
                                break;
                            default:
                                break;
                        }
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
                    if (step_no % stride_size != 0)
                    {
                        to_write = false;
                    }
                    if (to_write)
                    {
                        if (write_merged)
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

                if (write_parsed_data)
                {
                    fout_data << format("%1$-10.2f %2$-7d %3$-8d %4$14.4f %5$-20.12f %6$-20.12f %7$-20.12f") % step_time % step_no %
                    no_subsystems % step_temperature % step_potEne % step_kinEne % step_mdEne << endl; 
                }
            }
            else
            {
                buffer << str << endl;
            }
        }
        if (write_merged)
            fout << buffer.str() << endl;
        
        
        if (verbose)
        {
            cerr << min_step << ", " << max_step << endl;
        }
            
        last_step_nos.push_back(max_step);
    }
    return last_step_nos;
}

void merge_datafile(const std::vector<std::string> &folders, const vector<int>& last_step_nos,
    const std::string& filename, const std::string& merged_filename, int stride_size, bool verbose)
{
    string str;
    if (verbose)
        cout << "Combining */" << filename << " to " << merged_filename  << endl;

    ofstream fout( (merged_filename).c_str() );
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
            if (step_no % stride_size != 0)
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


void merge_mullfile(const std::vector<std::string> &folders, const vector<int>& last_step_nos,
    const std::string& filename, const std::string& merged_filename, int stride_size, bool convert_to_nac, bool verbose)
{
    string str;
    if (verbose)
        cout << "Combining */" << filename << " to " << merged_filename  << endl;

    ofstream fout( (merged_filename).c_str() );
    for (size_t i=0; i<folders.size(); ++i)
    {
        string filein = (folders[i] + "/" + filename );
        if (verbose)
            cout <<  "  Loading " << filein << endl;

        ifstream fin( filein.c_str() );

        bool first_line = true;

        int nat = 0;
        int norb = 0;
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
                norb = lexical_cast<int>(SplitVec[0]) ;  
                nat = lexical_cast<int>(SplitVec[1]) ;               
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
            if (step_no % stride_size != 0)
            {
                to_write = false;
            }
            if (to_write)
            {
                if (convert_to_nac)
                {
                    fout << nat << endl;
                }
                else
                {
                    fout << norb << endl;
                }               
                fout << title_line << endl;
            }
            if (convert_to_nac)
            {
                if (to_write)
                {
                    int cur_atom = 1;
                    string cur_sym;
                    double cur_charge = 0.0;
                    int ind_atom; 
                    string sym, orb;
                    double charge;
                    for (int j=0; j<norb;++j)
                    {
                        fin >> ind_atom >> sym >> orb >> charge;
                        if (ind_atom == cur_atom)
                        {
                            cur_charge += charge;
                            cur_sym = sym;
                        }
                        else
                        {
                            fout << format("%1$-4s %2$14.8f") % cur_sym % cur_charge << endl;
                            cur_charge = charge;
                            cur_atom = ind_atom;
                            cur_sym = sym;
                        }
                    }
                    fout << format("%1$-4s %2$14.8f") % cur_sym % cur_charge << endl;
                    getline(fin, str);
                }
                else
                {
                    for (int j=0; j<norb;++j)
                    {
                        getline(fin, str);
                    }
                }
            }
            else
            {
                for (int j=0; j<norb;++j)
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

    m.def("merge_mullfile", &merge_mullfile, R"pbdoc(
        Merge a specific mulliken file of DCDFTBMD.

        
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
