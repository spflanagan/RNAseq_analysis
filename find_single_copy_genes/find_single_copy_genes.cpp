//Author: Sarah P. Flanagan
//Date: 31 March 2015
//Purpose: Read in either a list of .fasta files or .paml files and print out a list of those files that have single copy orthologs. 
//Usage:
//This program can be run in either intreactive mode or on the command line 
//Find orthologs with single copy genes (./find_single_copy_genes):\n";
//-i: Input list (include path)
//-o: Output list (include path)
//-d: Directory containing sequence files listed in the input file
//-p: Add this flag if the sequences are in paml format instead of fasta. Default is fasta.
//-h: Print help message
//no arguments: interactive mode

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cstring>

using namespace std;

bool FileTest(ifstream& file, string filename)
{
	cout << filename;
	if (file.is_open())
		cout << " open\n";
	else
	{
		while (!file.is_open())
		{
			cout << " not open. Please re-enter filename: ";
			getline(cin, filename, '\n');
			file.open(filename);
		}
	}
	return true;
}

istream& universal_getline(istream& is, string& t)
{
	//this code is adapted from a post on stackoverflow:
	// http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
	//written by user763305
	t.clear();
	istream::sentry se(is, true);
	streambuf* sb = is.rdbuf();//sets pointer to stream buffer object

	for (;;)
	{
		int c = sb->sbumpc();//get current character and advance to the next position
		switch (c)//tests for equality against a list of variables (like multiple if statements)
		{
		case '\n'://if the next character is '\n', return the line
			return is;
		case '\r'://if the character is '\n', see if the next one is '\n' or return the line
			if (sb->sgetc() == '\n')
				sb->sbumpc();
			return is;
		case EOF:// Also handle the case when the last line has no line ending
			if (t.empty())//if there's nothing there, set it to be the end of file
				is.setstate(ios::eofbit);//set it to be the end of the file and return it
			return is;
		default://if none of the above, continue on.
			t += (char)c;
		}
	}

}

class fasta_record
{
public:
	string seq_id;
	string sequence;

	fasta_record()
	{
		seq_id = "";
		sequence = "";
	}
};

int read_fasta_record(ifstream &fasta_file, vector <fasta_record> &record)
{
	int count = 0;
	string line;
	while (!fasta_file.eof())
	{
		if (!fasta_file.eof())
		{
			universal_getline(fasta_file, line);
			if (line.substr(0, 1) == ">")
			{
				//it's the name
				record.push_back(fasta_record());
				count++;
				record[count - 1].seq_id = line.substr(1, line.size());
			}
			//if (line.substr(0, 1) == "A" || line.substr(0, 1) == "T" || line.substr(0, 1) == "C" || line.substr(0, 1) == "G" || line.substr(0, 1) == "N")
			else
			{//it's the sequence
				if (line.substr(0, 1) != "\n")
					record[count - 1].sequence.append(line);
			}
		}
	}
	cout << "Fasta file read!\n";
	return count;
}

int read_paml_record(ifstream &paml_file, vector <fasta_record> &record)
{
	int count = 0;
	string line, seq_name, seq;
	while (!paml_file.eof())
	{
		if (!paml_file.eof())
		{
			if (count == 0)
			{
				universal_getline(paml_file, line);
				count++;
			}
			if (count != 0)
			{
				universal_getline(paml_file, line);
				if (line.substr(0, 1) != "\n");
				{
					istringstream ss(line);
					ss >> seq_name;
					record.push_back(fasta_record());
					record[count - 1].seq_id = seq_name;
					ss >> seq;
					record[count - 1].sequence.append(seq);
					count++;
				}
			}	
		}
	}
	cout << "paml file read!\n";
	return count - 1;
}

int main(int argc, char* argv[])
{
	int i, end, num_sequences, count;
	string input_list_name, seq_file_name, output_list_name, seq_dir_name;
	ifstream input_list, seq_file;
	ofstream output_list;
	string line, last_species;
	bool fasta = true;
	bool single_copy = true;
	bool interactivemode = false;
	string tempstring1, tempstring2, query;
	
	seq_dir_name = "";//"E://ubuntushare//orthomcl_dir//paml//";
	input_list_name = "default";// "E://ubuntushare//orthomcl_dir//paml//test_seqs.txt";
	output_list_name = "single_copy.txt";// E://ubuntushare//orthomcl_dir//single_copy.txt";

	if (argc == 1)
	{
		cout << "\n(I)nteractive or (H)elp?\n";
		cin >> query;
		if (query == "H" || query == "h")
		{
			cout << "\nFind orthologs with single copy genes (./find_single_copy_genes):\n";
			cout << "-i:\tInput list (include path)\n";
			cout << "-o:\tOutput list (include path)\n";
			cout << "-d:\tProvide directory for sequences files in input list\n";
			cout << "-p:\tAdd this flag if the input is paml format instead of fasta format\n";
			cout << "-h:\tPrint help message\n";
			cout << "no arguments:\tinteractive mode\n";
			return 0;
		}
	}

	if (argc > 1)
	{
		tempstring1 = argv[1];
		if (tempstring1 == "-h")
		{
			cout << "\nFind orthologs with single copy genes (./find_single_copy_genes):\n";
			cout << "-i:\tInput list (include path)\n";
			cout << "-o:\tOutput list (include path)\n";
			cout << "-d:\tProvide directory for sequences files in input list\n";
			cout << "-p:\tAdd this flag plus a zero (-p 0) if the input is paml format instead of fasta format\n";
			cout << "-h:\tPrint help message\n";
			cout << "no arguments:\tinteractive mode\n";
			return 0;
		}
	}

	for (i = 1; i < argc - 1; i++)
	{
		interactivemode = false;
		tempstring1 = argv[i];
		tempstring2 = argv[i + 1];
		if (tempstring1 == "-i")
			input_list_name = tempstring2;
		if (tempstring1 == "-o")
			output_list_name = tempstring2;
		if (tempstring1 == "-d")
			seq_dir_name = tempstring2;
		if (tempstring1 == "-p")
			fasta = false;
	}

	
	if (interactivemode)
	{
		cout << "Provide input list file name (with the path):\n";
		cin >> input_list_name;
		cout << "Default output file name is single_copy.txt in the current directory.\nIs that acceptable? Y or N\n";
		cin >> tempstring1;
		if (tempstring1 == "N" | tempstring1 == "n")
		{
			cout << "Enter desired input directory:\t";
			cin >> output_list_name;
		}
		cout << "Provide directory for sequence files listed in the input file:\n";
		cin >> seq_dir_name;
		cout << "Default format of sequence files listed in input list is fasta format. Are your files paml format instead? Y or N\n";
		cin >> tempstring1;
		if (tempstring1 == "Y" | tempstring1 == "y")
		{
			fasta = false;
		}
	}
	if (input_list_name == "default")
	{
		cout << "\nProvide input list name\n";
		cin >> input_list_name;
		interactivemode = true;
	}
	if (seq_dir_name == "")
	{
		cout << "\nProvide directory for sequence files listed in input list file\n";
		cin >> seq_dir_name;
		interactivemode = true;
	}


	cout << "\n\nInput list file name:\t" << input_list_name;
	cout << "\nOutput list file name:\t" << output_list_name;
	cout << "\nDirectory containing sequence files:\t" << seq_dir_name;
	if (fasta)
		cout << "\nThe files listed in the input are fasta files.\n";
	else
		cout << "\nThe files listed in the input are paml files.\n";

	if (interactivemode)
	{
		cout << "\n\nProceed? (y to proceed)\n";
		cin >> query;

		if (query != "y" && query != "Y")
		{
			cout << "\n\nEnter an integer to exit!!\n";
			cin >> i;
			return 0;
		}
	}

	cout << "\n\nProceeding...\n";

	input_list.open(input_list_name);
	FileTest(input_list, input_list_name);
	output_list.open(output_list_name);
	count = 0;
	while (universal_getline(input_list, line))
	{
		if (!input_list.eof())
		{
			if (line != "")
			{
				seq_file_name = seq_dir_name + line;
				seq_file.open(seq_file_name);
				FileTest(seq_file, seq_file_name);
				vector<fasta_record> sequences;
				if (fasta)
					num_sequences = read_fasta_record(seq_file, sequences);
				if (!fasta)//then it's paml
					num_sequences = read_paml_record(seq_file, sequences);
				seq_file.close();
				single_copy = true;
				last_species = sequences[0].seq_id.substr(0, 4);
				for (i = 1; i < num_sequences; i++)
				{
					if (sequences[i].seq_id.substr(0, 4) == last_species)
						single_copy = false;
				}
				if (single_copy == true)
				{
					if (count == 0)
						output_list << line;
					else
						output_list << '\n' << line;
					count++;
				}
			}
		}
	}

	input_list.close();
	output_list.close();

	if (interactivemode)
	{
		cout << "Done! Enter integer to quit:\n";
		cin >> end;
		return 0;
	}
	else
	{
		cout << "Done!\n";
		return 0;
	}
}