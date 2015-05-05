//Author: Sarah P. Flanagan
//Date: 30 March 2015
//Purpose: Take a fasta file and convert it to a paml file. 
//Usage:
//This program can be run in either intreactive mode or on the command line 
//Convert fasta file format to input format for PAML (./fasta_to_paml):\n";
//-b: Base file name (e.g. ovary1000)
//-d: Direcotry name (include path)
//-o: Output Directory name and path (if not given, assumed same as input directory)
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

int main(int argc, char* argv[])
{
	int i, end, num_seqs, count, seq_length, spp_abbr_length;
	string directory_name, base_file_name, fasta_file_name, out_file_name, out_dir_name;
	ifstream fasta_file;
	ofstream out_file;
	vector <fasta_record> protein;
	string seq_name, line, line2;
	bool interactivemode = false;
	string tempstring1, tempstring2, query;
	stringstream ss;

	directory_name = "";
	base_file_name = "default";
	out_dir_name = directory_name;

	if (argc == 1)
	{
		cout << "\n(I)nteractive or (H)elp?\n";
		cin >> query;
		if (query == "H" || query == "h")
		{
			cout << "\nConvert fasta file format to input format for PAML (./fasta_to_paml):\n";
			cout << "-b:\tBase file name (e.g. ovary1000)\n";
			cout << "-d:\tInput Directory name (include path). If not provided, the files will write to the directory where this program is saved.\n";
			cout << "-o:\tOutput Directory name and path (if not given, assumed same as input directory)\n";
			cout << "-h:\tPrint help message\n";
			cout << "no arguments:\tinteractive mode\n";
			return 0;
		}
		if (query == "I" || query == "i")
		{
			interactivemode = true;
		}
	}

	if (argc > 1)
	{
		tempstring1 = argv[1];
		if (tempstring1 == "-h")
		{
			cout << "\nConvert fasta file format to input format for PAML (./fasta_to_paml):\n";
			cout << "-b:\tBase file name (e.g. ovary1000)\n";
			cout << "-d:\tInput Direcotry name (include path). If not provided, the files will write to the directory where this program is saved.\n";
			cout << "-o:\tOutput Directory name and path (if not given, assumed same as input directory)\n";
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
		if (tempstring1 == "-d")
			directory_name = tempstring2;
		if (tempstring1 == "-b")
			base_file_name = tempstring2;
		if (tempstring1 == "-o")
			out_dir_name = tempstring2;
	}

	
	if (interactivemode)
	{
		cout << "Input base file name:\n";
		cin >> base_file_name;
		cout << "Do you want to provide an input directory? Y or N\n";
		cin >> tempstring1;
		if (tempstring1 == "Y" | tempstring1 == "y")
		{
			cout << "Enter desired input directory:\t";
			cin >> directory_name;
		}
		cout << "Do you want to provide an output directory? Y or N\n";
		cin >> tempstring1;
		if (tempstring1 == "Y" | tempstring1 == "y")
		{
			cout << "Enter desired output directory:\t";
			cin >> out_dir_name;
		}
	}

	if (base_file_name == "default")
	{
		cout << "\nBase file name (e.g. ovary1000)\n";
		cin >> base_file_name;
		interactivemode = true;
	}

	cout << "\n\nBase file name:\t" << base_file_name;
	cout << "\nInput Directory name:\t" << directory_name;
	cout << "\nOutput Directory name:\t" << out_dir_name;

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

	ss.str(string());
	if (directory_name == "")
		ss << base_file_name << ".fasta";
	else
		ss << directory_name << base_file_name << ".fasta";
	fasta_file_name = ss.str();
	fasta_file.open(fasta_file_name);
	FileTest(fasta_file, fasta_file_name);
	num_seqs = read_fasta_record(fasta_file, protein);
	fasta_file.close();
	cout << num_seqs << " sequences found\n";
	seq_length = protein[0].sequence.length();
	for (i = 0; i < num_seqs; i++)
	{
		if (seq_length != protein[i].sequence.length())//make sure all the sequences are the same length
		{
			cout << fasta_file_name << ":\tNot all sequences are the same length! Cannot write in PAML format. Exiting program.\n";
			if (interactivemode)
			{
				cout << "Input integer to quit:\n";
				cin >> end;
				return 0;
			}
			else
			{
				return 0;
			}
		}
		for (size_t t = 0; t < protein[i].sequence.size()-3; t=t+3)
		{
			if (protein[i].sequence.substr(t, t + 3) == "TAG" || protein[i].sequence.substr(t, t + 3) == "tag" ||
				protein[i].sequence.substr(t, t + 3) == "TAA" || protein[i].sequence.substr(t, t + 3) == "taa" ||
				protein[i].sequence.substr(t, t + 3) == "TGA" ||protein[i].sequence.substr(t, t + 3) == "tga")//and remove stop codons
			{
				protein[i].sequence.replace(t, 3, "???");
				cout << protein[i].seq_id << " had stop codon.\n";
			}
		}
		
	}
	if (protein[0].seq_id.substr(3, 1) == "|")
		spp_abbr_length = 3;
	else
		spp_abbr_length = 4;
	
	ss.str(string());
	if (out_dir_name == "")
		ss << base_file_name << ".paml";
	else
		ss << out_dir_name << base_file_name << ".paml";
	out_file_name = ss.str();
	cout << "Writing to file " << out_file_name << '\n';
	ss.str(string());
	out_file.open(out_file_name);
	count = 0;
	for (i = 0; i < num_seqs; i++)
	{
		if (count == 0)
			out_file << num_seqs << " " << protein[i].sequence.length();
		out_file << '\n' << protein[i].seq_id.substr(0, spp_abbr_length) << "  " << protein[i].sequence;
		count++;
	}
	out_file.close();	

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