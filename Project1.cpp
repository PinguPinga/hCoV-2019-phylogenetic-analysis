#define _CRT_SECURE_NO_WARNINGS
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#define P_len 20
#include <map>

int fail_thres, mut_thres;
bool hq_nt, kind_taiwan;

inline int nt_num(char c) {
	return c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : 3;
}

inline int mut_num(int loc, char ori, char mut) {
	return ((((loc << 2) + nt_num(ori) ) << 2) + nt_num(mut));
}

inline char num_nt(int n) {
	return n == 0 ? 'A' : n == 1 ? 'C' : n == 2 ? 'G' : 'T';
}

inline std::string num_mut(int num) {
	int loc = num >> 4;
	char ori = num_nt((num >> 2) & 3), mut = num_nt(num & 3);
	std::stringstream ss;
	ss << ori << loc << mut;
	return ss.str();
}

struct pt {
	// >hCoV-19/USA/WA-S88/2020|EPI_ISL_417141|2020-03-01
	std::string seq; // sequence
	std::string acc; // accession id
	std::string virus_name, host_species, area_name = "", lab_id; // info
	int y = 0, m = 0, d = 0; // year, month, day
	pt() {}
	pt(std::string inf_in, std::string seq_in) {
		seq = seq_in;
		int idx1 = inf_in.find('|'), idx2 = inf_in.find('|', idx1 + 1);
		acc = inf_in.substr(idx1 + 1, idx2 - idx1 - 1);
		for(int i = 0; i < idx1; ++i) {
			if(inf_in[i] == ' ') // replace space with underline
				inf_in[i] = '_';
			else if(inf_in[i] == '/') // replace '/' with space
				inf_in[i] = ' ';
		}
		std::stringstream ss;
		ss << inf_in.substr(1, idx1 - 1);
		std::string tmp;
		ss >> virus_name >> tmp;
		if('a' <= tmp[0] && tmp[0] <= 'z') {
			host_species = tmp;
			ss >> tmp;
		}
		else
			host_species = "human";
		for(int i = 0; i < tmp.length(); ++i)
			if(isalpha(tmp[i]))
				area_name += tolower(tmp[i]);
		ss >> lab_id >> tmp;
		ss.clear();
		for(int i = idx2 + 1; i < inf_in.length(); ++i) {
			if(inf_in[i] == '-') // replace '-' with space
				inf_in[i] = ' ';
		}
		ss << inf_in.substr(idx2 + 1, 10);
		y = m = d = 0;
		ss >> y >> m >> d; // 2020 03 01
		ss.clear();
	}
};

struct orf {
	std::string name;
	int bg, ed;
	orf() {}
	orf(std::string name_in, int bg_in, int ed_in) {
		name = name_in;
		bg = bg_in;
		ed = ed_in;
	}
};

bool input_data(pt &ref, std::vector<orf> &orfs, std::vector<pt> &dat) {
	std::string buf; // buffer
	// read reference genomic data
	std::fstream ref_file; // reference (NC_045512.2)
	ref_file.open("SARS-CoV-2 (NC_045512.2).fasta", std::fstream::in);
	if(!ref_file.is_open()) {
		std::cout << "reference genome data reading error" << std::endl;
		return 0;
	}
	getline(ref_file, buf);
	ref = pt(buf, "");
	while(!ref_file.eof()) {
		getline(ref_file, buf);
		for(int i = 0; i < buf.length(); ++i)
			buf[i] = toupper(buf[i]);
		ref.seq += buf;
	}
	ref_file.close();
	// read ORF data
	std::fstream orf_file; // ORF
	orf_file.open("Project1_ORF_file.txt", std::fstream::in);
	if(!orf_file.is_open()) {
		std::cout << "reference genome data reading error" << std::endl;
		return 0;
	}
	while(!orf_file.eof()) {
		orf orf_tmp;
		orf_file >> orf_tmp.name >> orf_tmp.bg >> orf_tmp.ed;
		orfs.emplace_back(orf_tmp);
	}
	orf_file.close();
	// read GISAID fasta genome data
	std::fstream fst_file; // GISAID fasta file
	fst_file.open("gisaid_cov2020_sequences.fasta", std::fstream::in);
	if(!fst_file.is_open()) {
		std::cout << "GISAID genome database reading error" << std::endl;
		return 0;
	}
	while(!fst_file.eof()) { // while(!fst_file.eof() && dat.size() < 51) { // while(!fst_file.eof()) { // 
		getline(fst_file, buf);
		if(buf[0] == '>')
			dat.emplace_back(pt(buf, ""));
		else {
			for(int i = 0; i < buf.length(); ++i)
			buf[i] = toupper(buf[i]);
			dat[dat.size() - 1].seq += buf;
		}
	}
	fst_file.close();
	// dat.pop_back(); // for debug
	return 1;
}

void kmp_initialize_table(std::string &pat, std::vector<int> &kmp_tab) {
	int len = pat.length();
	kmp_tab[0] = - 1;
	int k = 0;
	for(int i = 1; i < len; ++i) {
		k = kmp_tab[i - 1];
		while(pat[i] != pat[k + 1]) {
			if(k == - 1)
				break;
			else
				k = kmp_tab[k];
		}
		if(pat[i] == pat[k + 1])
			kmp_tab[i] = k + 1;
		else
			kmp_tab[i] = - 1;
	}
}

int kmp_match(std::string &str, int bg, int ed, std::string &pat, std::vector<int> &tab) {
	for(int i = bg, j = 0; j < pat.length() && i < ed; ++i) {
		while(1) { // while(j < pat.length()) {
			if(str[i] == pat[j]) {
				++j;
				if(j == pat.length()) // str[i + 1 - pat.length() ... i] == pat, pattern found
					return i + 1 - pat.length();
				break;
			}
			else if(j == 0)
				break;
			else
				j = tab[j - 1] + 1;
		}
	}
	return -1;
}

bool find_mutation(pt &ref, std::vector<orf> &orfs, pt &cmp, std::vector<int> &loc, std::vector<std::string> &pri, std::vector<std::vector<int>> &kmp_tab, std::vector<int> &pms) {
	int fail = 0;
	std::vector<int> sft(orfs.size(), INT_MAX), lft(orfs.size(), INT_MAX), rgt(orfs.size(), INT_MAX);
	for(int i = 0; i < orfs.size(); ++i) {
		int dif = abs((int)(cmp.seq.length() - ref.seq.length()));
		int kmp_loc = kmp_match(cmp.seq, std::max(0, orfs[i].bg - 1 - dif), std::min((int)cmp.seq.length(), orfs[i].ed + dif), pri[i], kmp_tab[i]);
		if(kmp_loc != -1)
			sft[i] = kmp_loc - loc[i];
	}
	for(int i = 1; i < orfs.size(); ++i) {
		if(sft[i - 1] != INT_MAX)
			lft[i] = sft[i - 1];
		else if(lft[i - 1] != INT_MAX)
			lft[i] = lft[i - 1];
	}
	for(int i = orfs.size() - 2; i >= 0; --i) {
		if(sft[i + 1] != INT_MAX)
			rgt[i] = sft[i + 1];
		else if(lft[i + 1] != INT_MAX)
			rgt[i] = rgt[i + 1];
	}
	for(int i = 0; i < orfs.size(); ++i) {
		if(sft[i] == INT_MAX) {
			if(lft[i] != INT_MAX && (lft[i] == rgt[i] || i == orfs.size() - 1))
				sft[i] = lft[i];
			else if(rgt[i] != INT_MAX && i == 0)
				sft[i] = rgt[i];
			else {
				// std::cout << orfs[i].name << " not found, probably due to primer / sequence error. " << std::endl;
				++fail;
				continue;
			}
		}
		for(int j = /* orfs[i].bg == 266 ? 241 - 1 :*/ orfs[i].bg - 1; j < orfs[i].ed && j < ref.seq.length() && j + sft[i] < cmp.seq.length(); ++j) {
			if(cmp.seq[j + sft[i]] != 'A' && cmp.seq[j + sft[i]] != 'C' && cmp.seq[j + sft[i]] != 'G' && cmp.seq[j + sft[i]] != 'T') {
				if(hq_nt && !(kind_taiwan && cmp.area_name == "taiwan"))
					return 0;
				continue;
			}
			if(ref.seq[j] != cmp.seq[j + sft[i]]) { // no mutations important before 289
				// std::cout << orfs[i].name << " mutation " << ref.seq[j] << j+1 << cmp.seq[j + sft[i]] << " found" << std::endl;
				pms.emplace_back(mut_num(j + 1, ref.seq[j], cmp.seq[j + sft[i]]));
			}
		}
	}
	return (fail < fail_thres && pms.size() < mut_thres);
}

bool input_setting() {
	std::fstream setting_file;
	setting_file.open("Project1_setting_file.txt", std::fstream::in);
	if(!setting_file.is_open()) {
		std::cout << "no setting file" << std::endl;
		return 0;
	}
	std::string buf;
	std::stringstream ss;
	getline(setting_file, buf);
	ss.str(buf);
	ss >> fail_thres; // fail if >= N ORFs not detected
	getline(setting_file, buf);
	ss.str(buf);
	ss >> mut_thres; // fail if >= N mutations detected
	getline(setting_file, buf);
	ss.str(buf);
	ss >> hq_nt; // reject nucleotides except ACGT
	getline(setting_file, buf);
	ss.str(buf);
	ss >> kind_taiwan; // reject nucleotides except ACGT
	setting_file.close();
	return 1;
}

int main()
{
	if(input_setting() == 0) { // input
		system("pause");
		return 0;
	}
	pt ref; // reference genome
	std::vector<orf> orfs; // orf name and location according to reference genome
	std::vector<pt> dat; // GISAID genome data
	if(input_data(ref, orfs, dat) == 0) { // input
		system("pause");
		return 0;
	}
	// find_primer(dat, pri, loc); // find primer
	std::vector<std::string> pri(orfs.size()); // primer to find ORF
	std::vector<int> loc(orfs.size()); // primer location in sequence of reference genome
	for(int i = 0; i < orfs.size(); ++i) {
		loc[i] = orfs[i].bg - 1;
		pri[i] = ref.seq.substr(loc[i], P_len);
	}
	// initialize KMP table for pattern matching
	std::vector<std::vector<int>> kmp_tab(orfs.size(), std::vector<int>(P_len, 0)); // kmp table
	for(int i = 0; i < orfs.size(); ++i)
		kmp_initialize_table(pri[i], kmp_tab[i]);
	// output mutation start
	std::fstream mut_file, out_file;
	mut_file.open("hCoV-2019 point mutation data.txt", std::fstream::out);
	out_file.open("Project2_input_file.txt", std::fstream::out); // input file for Project 2
	std::vector<int> pms;
	for(int i = 0; i < dat.size(); ++i) {
		mut_file << dat[i].acc << '	' << dat[i].area_name << '	' << dat[i].y << '-' << std::setfill('0') << std::setw(2) << dat[i].m << '-' << std::setfill('0') << std::setw(2) << dat[i].d << std::endl;
		if(find_mutation(ref, orfs, dat[i], loc, pri, kmp_tab, pms)) {
			out_file << dat[i].acc << '	' << dat[i].area_name << '	' << dat[i].y << '	' << dat[i].m << '	' << dat[i].d << std::endl;
			out_file << pms.size();
			for(int j = 0; j < pms.size(); ++j) {
				mut_file << num_mut(pms[j]);
				out_file << '	' << pms[j];
				if(j != pms.size() - 1)
					mut_file << '	';
			}
			mut_file << std::endl;
			out_file << std::endl << std::endl;
		}
		else
			mut_file << "WARNING: poor sequence quality" << std::endl;
		mut_file << std::endl;
		pms.clear();
	}
	mut_file.close();
	out_file.close();
	return 0;
}