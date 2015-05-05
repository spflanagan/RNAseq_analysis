//Author: Sarah P. Flanagan
//Date: 27 March 2015
//Purpose: Take the groups.txt file output by orthomcl and pull out all of the sequences 
//Usage:
//(I)nteractive or (H)elp?
//Convert OrthoMCL 'groups.txt' output to fasta (./om_groups_to_fasta)
//Converts orthomcl output to fasta format
//-g: Groups file, including path
//-f: Fasta file name, including path. Can provide multiples, but must be separated by commas
//-n: If this is the original nucleotide file without corrected headers, follow -n flag with species code. If it has the appropriate header, input false
//-o: output directory name
//-h: Output the usage message
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

int read_fasta_record(ifstream &fasta_file, vector <fasta_record> &record, string spp)
{
	int count = 0;
	string line, new_name;
	spp = spp + "|";
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
				istringstream seq(line.substr(1, line.size()));
				seq >> new_name;
				record[count - 1].seq_id = spp + new_name;
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
	int end, count, i, num_protein_seqs, group_count;
	size_t t;
	string groups_file_name, out_file_dir, out_file_name, fasta_name, abbr;
	vector<string> fasta_file_names, spp_abbr;
	ifstream groups_file, fasta_file;
	ofstream fasta_names;
	vector <fasta_record> protein;
	string seq_name, line, line2;
	bool wrong_header = false;
	bool interactivemode = true;
	string tempstring1, tempstring2, query;


	groups_file_name = "groups.txt";//"E:\\ubuntushare\\orthomcl_dir\\groups.txt";
	out_file_dir = "";//"E:\\ubuntushare\\orthomcl_dir\\groups\\";

	if (argc == 1)
	{
		cout << "\n(I)nteractive or (H)elp?\n";
		cin >> query;
		if (query == "H" || query == "h")
		{
			cout << "\nConvert OrthoMCL 'groups.txt' output to fasta (./om_groups_to_fasta):\n";
			cout << "Converts orthomcl output to fasta format\n";
			cout << "-g:\tGroups file, including path\n";
			cout << "-f:\tFasta file name, including path. Can provide multiples, but must be separated by commas\n";
			cout << "-n:\tIf this is the original nucleotide file without corrected headers, follow -n flag with species code. If it has the appropriate header, input false\n\t  Separate multiple entries with a comma";
			cout << "-o:\toutput directory name\n";
			cout << "-h:\tPrint this message\n";
			cout << "no arguments:\tinteractive mode\n";
			return 0;
		}
	}

	if (argc > 1)
	{
		tempstring1 = argv[1];
		if (tempstring1 == "-h")
		{
			cout << "\nConvert OrthoMCL 'groups.txt' output to fasta (./om_groups_to_fasta):\n";
			cout << "Converts orthomcl output to fasta format\n";
			cout << "-g:\tGroups file, including path\n";
			cout << "-f:\tFasta file name(s), including path. Can provide multiples, but must be separated by commas\n";
			cout << "-n:\tIf this is the original nucleotide file without corrected headers, follow -n flag with species code. If it has the appropriate header, input false.\n\t  Separate multiple entries with a comma";
			cout << "-o:\toutput directory name\n";
			cout << "no arguments:\tinteractive mode\n";
			return 0;
		}
	}

	for (i = 1; i < argc - 1; i++)
	{
		interactivemode = false;
		tempstring1 = argv[i];
		tempstring2 = argv[i + 1];
		if (tempstring1 == "-g")
			groups_file_name = tempstring2;
		if (tempstring1 == "-f")
			fasta_file_names.push_back(tempstring2);
		if (tempstring1 == "-n")
			abbr = tempstring2;
		if (abbr != "f" || abbr != "F" || abbr != "False" || abbr != "false" || abbr != "FALSE")
			spp_abbr.push_back(abbr);
		if (tempstring1 == "-o")
			out_file_dir = tempstring2;
	}

	if (interactivemode)
	{
		cout << "Default groups file name is " << groups_file_name << ". Is that acceptable? Y or N\n";
		cin >> tempstring1;
		if (tempstring1 == "N" || tempstring1 == "n")
		{
			cout << "Provide groups file name, including path\n";
			cin >> groups_file_name;
		}
		cout << "Provide Fasta file name(s). If you have multiples, separate with a comma:\n";
		cin >> fasta_name;
		istringstream names(fasta_name);
		while (getline(names, line, ','))
			fasta_file_names.push_back(line);
		cout << "Does your fasta file have the correct header? Y or N\n";
		cin >> tempstring1;
		if (tempstring1 == "N" || tempstring1 == "n")
		{
			cout << "Provide species abbreviation to convert header\n";
			cin >> abbr;
			istringstream names(abbr);
			while (getline(names, line, ','))
				spp_abbr.push_back(line);
			wrong_header = true;
		}
		cout << "Provide output directory:\n";
		cin >> out_file_dir;
	}

	if (groups_file_name == "groups.txt")
	{
		cout << "\nGroups file, including path\n";
		cin >> groups_file_name;
		interactivemode = true;
	}
	if (fasta_file_names.size() == 0)
	{
		cout << "\nProvide at least one fasta file name (separate multiples with a comma)\n";
		cin >> fasta_name;
		istringstream names(fasta_name);
		while (getline(names, line, ','))
			fasta_file_names.push_back(line);
		interactivemode = true;
	}
	if (out_file_dir == "")
	{
		cout << "\nInput output directory name\n\n";
		cin >> out_file_dir;
		interactivemode = true;
	}
	if (spp_abbr.size() == 0)
		wrong_header = false;
	else
		wrong_header = true;
	cout << "\n\nGroups file, including path:\t" << groups_file_name;
	for (t = 0; t < fasta_file_names.size(); t++)
	{
		cout << "\nFasta file, including path:\t" << fasta_file_names[t];
		if (!wrong_header)
			cout << "\nThis file has the correct header";
		else
			cout << "\nThis file has the incorrect header, which will be corrected using species abbreviation " << spp_abbr[t];
	}
	cout << "\nOutput directory name:\t" << out_file_dir;
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
	else
		cout << "\n\nProceeding...\n";

	num_protein_seqs = 0;
	for (t = 0; t < fasta_file_names.size(); t++)
	{
		fasta_file.open(fasta_file_names[t]);
		FileTest(fasta_file, fasta_file_names[t]);
		if (wrong_header)
			num_protein_seqs = num_protein_seqs + read_fasta_record(fasta_file, protein, spp_abbr[t]);
		else
			num_protein_seqs = num_protein_seqs + read_fasta_record(fasta_file, protein, "");
		fasta_file.close();
		cout << num_protein_seqs << " sequences found\n";
	}
	fasta_names.open("all_fasta_names.txt");
	for (t = 0; t < protein.size(); t++)
	{
		fasta_names << protein[t].seq_id << '\n';
		if (protein[t].sequence == "" && protein[t].seq_id == "")
		{
			protein[t].sequence.erase();
			protein[t].seq_id.erase();
			protein.erase(protein.begin() + t);
			num_protein_seqs--;
		}
	}
	fasta_names.close();

	groups_file.open(groups_file_name);
	FileTest(groups_file, groups_file_name);
	group_count = 0;
	while (!groups_file.eof())
	{
		universal_getline(groups_file, line);
		istringstream ss(line);
		ss >> seq_name; //first component is the sequence name
		if (seq_name.substr(seq_name.size() - 1, seq_name.size()) == ":")
			out_file_name = out_file_dir + seq_name.substr(0, seq_name.size() - 1);
		else
			out_file_name = out_file_dir + seq_name;
		out_file_name = out_file_name + ".fasta";
		vector<string> group_names;
		while (getline(ss, line2, ' '))//read in all of the group names
			group_names.push_back(line2);
		for (t = 0; t < group_names.size(); t++)
		{
			if (group_names[t] == "")
				group_names.erase(group_names.begin() + t);
		}
		ofstream out_file;
		out_file.open(out_file_name);
		count = 0;
		cout << out_file_name << ":\n";
		for (t = 0; t < group_names.size(); t++)
		{
			cout << group_names[t] << '\n';
			for (i = 0; i < num_protein_seqs; i++)
			{
				if (group_names[t].compare(protein[i].seq_id) == 0)//match group names to protein sequences
				{
					cout << "\tfound: " << protein[i].seq_id << '\n';
					if (count == 0)
						out_file << ">" << group_names[t] << '\n' << protein[i].sequence;
					else
						out_file << '\n' << ">" << group_names[t] << '\n' << protein[i].sequence;
					count++;
				}
			}			
		}
		out_file.close();
		group_count++;
		if (group_count % 1000 == 0)
			cout << group_count << " groups processed.\n";
	}
	groups_file.close();

	if (interactivemode)
	{
		cout << "Done! Input integer to quit:\n";
		cin >> end;
		return 0;
	}
	else
	{
		cout << "Done!\n";
		return 0;
	}
}