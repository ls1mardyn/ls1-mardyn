//
// Created by alex on 22.05.23.
//

#include "utils/OptionParser.h"
#include "utils/FileUtils.h"
#include "VTKFile.h"
#include "Analysis.h"

#include <fstream>

/**
 * @brief Initialize command line options.
 */
void initOptions(optparse::OptionParser& op) {
    op.usage("%prog [OPTIONS] <configfilename>\n\n"
              "Use option --help to display all available options.\n\n"
              "Analyses the provided vtk file for density by component id.\n");
    op.version("%prog 1.0");
    op.description("vtk-analysis (computes component densities of vtk file)");
    op.add_option("-n", "--samples-per-dim").dest("samples-per-dim").type("int").metavar("NUM").set_default(10).help("set how many samples are taken in every dimensions (default: %default)");
    op.add_option("-w", "--box-width").dest("box-width").type("float").metavar("NUM").set_default(10.0).help("width of the box for each sample (default: %default)");
    op.add_option("--bbox0").dest("bbox0").type("float").metavar("NUM").set_default(100.0).help("size of bounding box in dim 0 of simulation (default: %default)");
    op.add_option("--bbox1").dest("bbox1").type("float").metavar("NUM").set_default(100.0).help("size of bounding box in dim 1 of simulation (default: %default)");
    op.add_option("--bbox2").dest("bbox2").type("float").metavar("NUM").set_default(100.0).help("size of bounding box in dim 2 of simulation (default: %default)");
    op.add_option("-o", "--output").dest("output").type("string").metavar("PATH").set_default("vtk_analysis.txt").help("output file name (default: %default)");
    op.add_option("--check0").dest("check0").type("int").set_default(1).help("enables domain slicing in dimension 0 (default: %default)");
    op.add_option("--check1").dest("check1").type("int").set_default(0).help("enables domain slicing in dimension 1 (default: %default)");
    op.add_option("--check2").dest("check2").type("int").set_default(0).help("enables domain slicing in dimension 2 (default: %default)");
}

void writeOutput(const std::string& name, const std::vector<Analysis::Density>& data) {
    std::ofstream file;
    file.open (name);
    if(!file.is_open()) {
        std::cerr << "Unable to open output file!" << std::endl;
        exit(-1);
    }

    file << "###########################################\n";
    file << "                  pos x                    \n";
    file << "###########################################\n";
    for(auto& density : data) {
        file << density.position[0] << " ";
    }
    file << "\n";
    file << "###########################################\n";
    file << "                  pos y                    \n";
    file << "###########################################\n";
    for(auto& density : data) {
        file << density.position[1] << " ";
    }
    file << "\n";
    file << "###########################################\n";
    file << "                  pos z                    \n";
    file << "###########################################\n";
    for(auto& density : data) {
        file << density.position[2] << " ";
    }
    file << "\n";
    file << "###########################################\n";
    file << "                 density                   \n";
    file << "###########################################\n";
    for(auto& density : data) {
        for(auto val : density.densityValues)
            file << val << " ";
    }
    file << std::endl;
    file.flush();

    file.close();
}

int main(int argc, char** argv) {
    optparse::OptionParser op;
    initOptions(op);
    optparse::Values options = op.parse_args(argc, argv);
    std::vector<std::string> args = op.args();

    auto numArgs = args.size();
    if(numArgs != 1) {
        std::cerr << "Incorrect number of arguments provided." << std::endl;
        op.print_usage();
        exit(-1);
    }

    int nSamples = options.get("samples-per-dim");
    float boxWidth = options.get("box-width");
    std::array<float,3> bbox = {options.get("bbox0"), options.get("bbox1"), options.get("bbox2")};

    VTKFile vtk;
    std::string inputFile(args[0]);
    if( fileExists(inputFile.c_str()) ) {
        vtk.readXML(inputFile);
    } else {
        std::cerr << "Cannot open config file '" << inputFile << "'" << std::endl;
        exit(-2);
    }

    Analysis ana(vtk.getMolecules());
    std::array<int, 3> check { options.get("check0"), options.get("check1"), options.get("check2")};
    std::vector<Analysis::Density> densities = ana.computeDensities(nSamples, boxWidth, bbox, check);


    std::string fileName = std::string(options.get("output"));
    writeOutput(fileName, densities);


    return 0;
}
