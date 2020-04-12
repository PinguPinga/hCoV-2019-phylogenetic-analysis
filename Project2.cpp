#define _CRT_SECURE_NO_WARNINGS
#include <string>
#include <vector>
#include <iomanip>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#define H 196613

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

struct inf {
	std::string area; // sequence
	std::string acc; // accession id
	int y = 0, m = 0, d = 0; // year, month, day
	std::vector<int> mut;
};

int hash(std::vector<int> &vec) {
	int h = 0;
	for(auto i : vec)
		h = (h * 13 + i) % H;
	return h;
}

double muts_p(std::unordered_map<int, double> &mut_p, std::vector<int> &muts) {
	double p = 1;
	for(int i = 0; i < muts.size(); ++i)
		p *= mut_p[muts[i]];
	return p;
}

std::string muts_str(std::vector<int> &muts) {
	if(muts.empty())
		return "NIL";
	std::string str;
	for(int i = 0; i < muts.size(); ++i) {
		str += num_mut(muts[i]);
		if(i != muts.size() - 1)
			str += '	';
	}
	return str;
}

void input_pt_inf(std::vector<inf> &pt_inf) {
	std::fstream input_file;
	input_file.open("input_file.txt", std::fstream::in);
	if(!input_file.is_open())
		std::cout << "input data reading error" << std::endl;
	int n = 0, muts;
	while(!input_file.eof()) {
		pt_inf.emplace_back(inf());
		input_file >> pt_inf[n].acc >> pt_inf[n].area >> pt_inf[n].y >> pt_inf[n].m >> pt_inf[n].d;
		if(pt_inf[n].acc == "") {
			pt_inf.pop_back();
			break;
		}
		if(pt_inf[n].m == 0)
			pt_inf[n].m = 99;
		if(pt_inf[n].d == 0)
			pt_inf[n].d = 99;
		for(int i = 0; i < pt_inf[n].area.length(); ++i)
			pt_inf[n].area[i] = tolower(pt_inf[n].area[i]);
		input_file >> muts;
		pt_inf[n].mut.resize(muts);
		for(int i = 0; i < muts; ++i)
			input_file >> pt_inf[n].mut[i];
		++n;
		/* if(n == 50)
			break; */
	}
	input_file.close();
}

int find_nd(std::vector<std::vector<int>> &nd_muts, std::unordered_map<int, std::vector<int>> &hash_nd, std::vector<int> &mut) {
	int h = hash(mut);
	if(hash_nd.find(h) != hash_nd.end())
		for(int i = 0; i < hash_nd[h].size(); ++i) {
			if(nd_muts[hash_nd[h][i]].size() != mut.size())
				continue;
			int j;
			for(j = 0; j < mut.size(); ++j)
				if(nd_muts[hash_nd[h][i]][j] != mut[j])
					break;
			if(j == mut.size())
				return hash_nd[h][i];
		}
	return -1;
}

bool add_node(std::vector<std::vector<int>> &nd_muts, 
	std::unordered_map<int, std::vector<int>> &hash_nd, 
	std::unordered_map<int, std::vector<int>> &mut_nd, 
	std::vector<int> &nd_par, 
	std::vector<std::unordered_set<int>> &nd_chd, 
	std::vector<int> &nd_cnt, 
	std::vector<int> &muts) {
	int nd = find_nd(nd_muts, hash_nd, muts);
	if(nd != -1)
		return 0;
	nd = (int) nd_muts.size();
	hash_nd[hash(muts)].emplace_back(nd);
	for(int i = 0; i < muts.size(); ++i)
		mut_nd[muts[i]].emplace_back(nd);
	nd_muts.emplace_back(muts);
	nd_par.emplace_back(-1);
	nd_chd.emplace_back(std::unordered_set<int>());
	nd_cnt.emplace_back(0);
	return 1;
}

void generate_nodes(std::vector<inf> &pt_inf, 
	std::vector<std::vector<int>> &nd_pt, 
	std::vector<std::vector<int>> &nd_muts, 
	std::unordered_map<int, std::vector<int>> &hash_nd, 
	std::unordered_map<int, std::vector<int>> &mut_nd, 
	std::vector<int> &nd_par, 
	std::vector<std::unordered_set<int>> &nd_chd, 
	std::vector<int> &nd_cnt) {
	for(int i = 0; i < pt_inf.size(); ++i) {
		add_node(nd_muts, hash_nd, mut_nd, nd_par, nd_chd, nd_cnt, pt_inf[i].mut);
		int nd = find_nd(nd_muts, hash_nd, pt_inf[i].mut);
		nd_pt.resize(nd_muts.size());
		nd_pt[nd].emplace_back(i);
		++nd_cnt[nd];
	}
}

void del_node(std::vector<std::vector<int>> &nd_muts, 
	std::unordered_map<int, std::vector<int>> &hash_nd, 
	std::unordered_map<int, std::vector<int>> &mut_nd, 
	std::vector<int> &nd_par, 
	std::vector<std::unordered_set<int>> &nd_chd, 
	std::vector<int> &nd_cnt, 
	std::vector<bool> &prn, int nd) {
	if(prn[nd] == 1)
		return;
	prn[nd] = 1;
	// std::cout << "delete node #" << nd << std::endl;
	std::vector<int> tmp_muts = nd_muts[nd];
	for(auto j : nd_muts[nd]) {
		if(mut_nd[j].size() == 1)
			tmp_muts.erase(std::find(tmp_muts.begin(), tmp_muts.end(), j));
		// mut_nd[j].erase(std::find(mut_nd[j].begin(), mut_nd[j].end(), nd));
	}
	add_node(nd_muts, hash_nd, mut_nd, nd_par, nd_chd, nd_cnt, tmp_muts);
	int tmp_nd = find_nd(nd_muts, hash_nd, tmp_muts);
	// nd_par[nd] = tmp_nd;
	nd_chd[tmp_nd].insert(nd);
	if(nd_cnt[tmp_nd] > 0) {
		nd_cnt[tmp_nd] += nd_cnt[nd];
		nd_par[nd] = tmp_nd;
		// nd_cnt[nd] = 0;
	}
}

void prune_nodes(std::vector<std::vector<int>> &nd_muts, 
	std::unordered_map<int, std::vector<int>> &hash_nd, 
	std::unordered_map<int, std::vector<int>> &mut_nd, 
	std::vector<int> &nd_par, 
	std::vector<std::unordered_set<int>> &nd_chd, 
	std::vector<int> &nd_cnt, 
	std::vector<bool> &prn) {
	std::vector<int> mut_lst(mut_nd.size());
	prn.resize(nd_muts.size(), 0);
	for(std::unordered_map<int, std::vector<int>>::iterator i = mut_nd.begin(); i != mut_nd.end(); ++i)
		if(i->second.size() == 1) // del node with this unique mutation
			del_node(nd_muts, hash_nd, mut_nd, nd_par, nd_chd, nd_cnt, prn, i->second[0]);
}

void back_track(std::vector<std::vector<int>> &nd_muts, 
	std::unordered_map<int, std::vector<int>> &hash_nd, 
	std::unordered_map<int, std::vector<int>> &mut_nd, 
	std::vector<int> &nd_par, 
	std::vector<std::unordered_set<int>> &nd_chd, 
	std::vector<int> &nd_cnt, 
	std::vector<bool> &prn, 
	std::unordered_map<int, double> &mut_p) {
	std::vector<std::vector<int>> nd_q;
	for(int i = 0; i < nd_muts.size(); ++i) { // node 0 is root {no mutation}
		if(nd_par[i] == -1 && nd_muts[i].size() != 0 && (i >= prn.size() || prn[i] == 0)) {
			if(nd_q.size() < nd_muts[i].size() + 1)
				nd_q.resize(nd_muts[i].size() + 1);
			nd_q[nd_muts[i].size()].emplace_back(i);
		}
	}
	for(int i = (int) nd_q.size() - 1; i > 0; --i) { // layer of queue = mutations
		sort(nd_q[i].begin(), nd_q[i].end(), [&](int &a, int &b) {
			double a_p = muts_p(mut_p, nd_muts[a]), b_p = muts_p(mut_p, nd_muts[b]);
			if(a_p != b_p)
				return a_p > b_p;
			return a < b;
		} );
		for(int j = 0; j < nd_q[i].size(); ++j) {
			int nd = nd_q[i][j];
			// std::cout << "process node #" << nd << std::endl;
			for(std::unordered_set<int>::iterator k = nd_chd[nd].begin(); k != nd_chd[nd].end();) // refresh child list
				if(nd_par[*k] != -1 && nd_par[*k] != nd)
					nd_chd[nd].erase(k++);
				else
					++k;
			if(nd_chd[nd].empty() && nd_cnt[nd] == 0) // no need to trace back if no child and non-existent
				continue;
			if(nd_chd[nd].size() > 1 && nd_cnt[nd] == 0) { // create common ancestor of >=2 nodes
				for(auto k : nd_chd[nd]) {
					nd_par[k] = nd;
					nd_cnt[nd] += nd_cnt[k];
					// nd_cnt[k] = 0;
				}
			}
			int lnk_nd = nd_cnt[nd] > 0 ? nd : *nd_chd[nd].begin();
			std::vector<int> lst(nd_muts[nd].size());
			for(int k = 0; k < lst.size(); ++k)
				lst[k] = k;
			sort(lst.begin(), lst.end(), [&] (int &a, int &b) {
				if(mut_p[a] != mut_p[b])
					return mut_p[a] < mut_p[b];
				return a < b;
			} );
			int tmp_nd;
			std::vector<int> tmp_muts(nd_muts[nd].size() - 1);
			for(int k = 0; k < lst.size(); ++k) {
				for(int l = 0; l < lst[k]; ++l)
					tmp_muts[l] = nd_muts[nd][l];
				for(int l = lst[k] + 1; l < nd_muts[nd].size(); ++l)
					tmp_muts[l - 1] = nd_muts[nd][l];
				tmp_nd = find_nd(nd_muts, hash_nd, tmp_muts);
				if(tmp_nd != -1 && nd_cnt[tmp_nd] > 0) { //link lnk_nd to tmp_nd
					nd_chd[tmp_nd].insert(lnk_nd);
					nd_par[lnk_nd] = tmp_nd;
					nd_cnt[tmp_nd] += nd_cnt[lnk_nd];
					// nd_cnt[lnk_nd] = 0;
					break;
				}
			}
			if(nd_par[lnk_nd] != -1)
				continue;
			for(int k = 0; k + 1 < nd_muts[nd].size(); ++k) // copy nd_muts[1 .. last]
				tmp_muts[k] = nd_muts[nd][k + 1];
			bool add = add_node(nd_muts, hash_nd, mut_nd, nd_par, nd_chd, nd_cnt, tmp_muts);
			tmp_nd = find_nd(nd_muts, hash_nd, tmp_muts);
			if(add == 1) // added new node, which is currently non-existent
				nd_q[tmp_muts.size()].emplace_back(tmp_nd);
			nd_chd[tmp_nd].insert(lnk_nd);
			for(int k = 0; k + 1 < nd_muts[nd].size(); ++k) {
				tmp_muts[k] = nd_muts[nd][k];
				add = add_node(nd_muts, hash_nd, mut_nd, nd_par, nd_chd, nd_cnt, tmp_muts);
				tmp_nd = find_nd(nd_muts, hash_nd, tmp_muts);
				if(add == 1) // added new node, which is currently non-existent
					nd_q[tmp_muts.size()].emplace_back(tmp_nd);
				nd_chd[tmp_nd].insert(lnk_nd);
			}
		}
		nd_q[i].clear();
	}
}

void output_parent(std::vector<inf> &pt_inf, 
	std::vector<std::vector<int>> &nd_pt, 
	std::vector<std::vector<int>> &nd_muts, 
	std::vector<int> &nd_par, 
	std::vector<int> &nd_cnt) {
	std::fstream par_file;
	par_file.open("parent_file.txt", std::fstream::out);
	for(int i = 0; i < nd_muts.size(); ++i) {
		if(nd_cnt[i] == 0)
			continue;
		if(i < nd_pt.size() && nd_par[i] != -1 && nd_muts[i].size() - nd_muts[nd_par[i]].size() > 2) // nd_pt > 1 ? or subtree > 1?
			par_file << "WARNING: multiple (>=3) additional mutations" << std::endl;
		if(i >= nd_pt.size())
			par_file << "generated ";
		par_file << "node #" << i << std::endl;
		if(i < nd_pt.size())
			for(int j = 0; j < nd_pt[i].size(); ++j)
				par_file << pt_inf[nd_pt[i][j]].acc << ' ' << pt_inf[nd_pt[i][j]].area << ' '<< pt_inf[nd_pt[i][j]].y << '-'<< std::setfill('0') << std::setw(2) << pt_inf[nd_pt[i][j]].m << '-'<< std::setfill('0') << std::setw(2) << pt_inf[nd_pt[i][j]].d << std::endl;
		par_file << muts_str(nd_muts[i]) << std::endl;
		if(nd_par[i] >= nd_pt.size())
			par_file << "generated ";
		par_file << "parent #" << nd_par[i] << std::endl;
		if(nd_par[i] == -1)
			par_file << "NIL" << std::endl;
		else
			par_file << muts_str(nd_muts[nd_par[i]]) << std::endl;
		if(i != nd_muts.size() - 1)
			par_file << std::endl;
	}
	par_file.close();
}

void dfs_newick(std::fstream &out_file, 
	std::vector<inf> &pt_inf, 
	std::vector<std::vector<int>> &nd_pt, 
	std::vector<std::vector<int>> &nd_muts, 
	std::vector<std::unordered_set<int>> &nd_chd, 
	std::vector<int> &nd_cnt, 
	std::unordered_map<int, double> &mut_p, 
	int nd) {
	if(nd_chd[nd].empty()) {
		out_file << nd;
		return;
	}
	std::vector<int> lst;
	for(std::unordered_set<int>::iterator i = nd_chd[nd].begin(); i != nd_chd[nd].end(); ++i)
		lst.emplace_back(*i);
	sort(lst.begin(), lst.end(), [&] (int &a, int &b) {
			if(nd_cnt[a] != nd_cnt[b])
				return nd_cnt[a] > nd_cnt[b];
			double a_p = muts_p(mut_p, nd_muts[a]), b_p = muts_p(mut_p, nd_muts[b]);
			if(a_p != b_p)
				return a_p > b_p;
			return a < b;
	} );
	out_file << '(';
	for(int i = 0; i < lst.size(); ++i) {
		dfs_newick(out_file, pt_inf, nd_pt, nd_muts, nd_chd, nd_cnt, mut_p, lst[i]);
		if(i != lst.size() - 1)
			out_file << ',';
	}
	out_file << ')';
	out_file << nd;
}

void dfs_txt(std::fstream &out_file, 
	std::vector<inf> &pt_inf, 
	std::vector<std::vector<int>> &nd_pt, 
	std::vector<std::vector<int>> &nd_muts, 
	std::vector<int> &nd_par, 
	std::vector<std::unordered_set<int>> &nd_chd, 
	std::vector<int> &nd_cnt, 
	std::unordered_map<int, double> &mut_p, 
	int nd, std::string prefix) {
	// std::cout << "node #" << nd << "	"; // << ":	" << muts_str(nd_muts[nd]) << "		";
	out_file << "node #" << nd << "	"; // << ":	" << muts_str(nd_muts[nd]) << "		";
	if(nd < nd_pt.size()) {
		// std::cout << pt_inf[nd_pt[nd][0]].acc << ' ' << pt_inf[nd_pt[nd][0]].area << ' ' << pt_inf[nd_pt[nd][0]].y << '-'<< std::setfill('0') << std::setw(2) << pt_inf[nd_pt[nd][0]].m << '-' << std::setfill('0') << std::setw(2) << pt_inf[nd_pt[nd][0]].d;
		int t = 0;
		for(int i = 1; i < nd_pt[nd].size(); ++i)
			if(pt_inf[nd_pt[nd][t]].y * 12 + pt_inf[nd_pt[nd][t]].m * 31 + pt_inf[nd_pt[nd][t]].d > pt_inf[nd_pt[nd][i]].y * 12 + pt_inf[nd_pt[nd][i]].m * 31 + pt_inf[nd_pt[nd][i]].d)
				t = i;
		out_file << pt_inf[nd_pt[nd][t]].acc << ' ' << pt_inf[nd_pt[nd][t]].area << ' ' << pt_inf[nd_pt[nd][t]].y << '-'<< std::setfill('0') << std::setw(2) << pt_inf[nd_pt[nd][t]].m << '-' << std::setfill('0') << std::setw(2) << pt_inf[nd_pt[nd][t]].d;
		if(nd_pt[nd].size() > 1 || nd_cnt[nd] > nd_pt[nd].size() + nd_chd[nd].size()) {
			out_file << " (";
			if(nd_pt[nd].size() == 1)
				out_file << "1 entry"; 
			else
				out_file << nd_pt[nd].size() << " entries"; 
			if(nd_cnt[nd] > nd_pt[nd].size() + nd_chd[nd].size())
				out_file << ", " << nd_cnt[nd] << " entries under clade";
			out_file << ")";
		}
	}
	if(nd_par[nd] != -1 && nd_muts[nd].size() - nd_muts[nd_par[nd]].size() > 2)
		out_file << "	WARNING: multiple (>=3) additional mutations";
	// std::cout << std::endl;
	out_file << std::endl;
	if(nd_chd[nd].empty())
		return;
	std::vector<int> lst;
	for(std::unordered_set<int>::iterator i = nd_chd[nd].begin(); i != nd_chd[nd].end(); ++i)
		lst.emplace_back(*i);
	sort(lst.begin(), lst.end(), [&] (int &a, int &b) {
			if(nd_cnt[a] != nd_cnt[b])
				return nd_cnt[a] > nd_cnt[b];
			double a_p = muts_p(mut_p, nd_muts[a]), b_p = muts_p(mut_p, nd_muts[b]);
			if(a_p != b_p)
				return a_p > b_p;
			return a < b;
	} );
	for(int i = 0; i < lst.size(); ++i) {
		// std::cout << prefix;
		out_file << prefix;
		if(i != lst.size() - 1) {
			// std::cout << "¢u";
			out_file << "¢u";
			dfs_txt(out_file, pt_inf, nd_pt, nd_muts, nd_par, nd_chd, nd_cnt, mut_p, lst[i], prefix + "¢x");
		}
		else if(i == lst.size() - 1) {
			// std::cout << "¢|";
			out_file << "¢|";
			dfs_txt(out_file, pt_inf, nd_pt, nd_muts, nd_par, nd_chd, nd_cnt, mut_p, lst[i], prefix + "¡@");
		}
	}
}

int main()
{
	std::vector<inf> pt_inf;
	std::vector<int> nd_par, nd_cnt, empty_vector;
	std::vector<std::vector<int>> nd_muts, nd_pt;
	std::vector<std::unordered_set<int>> nd_chd;
	std::unordered_map<int, std::vector<int>> hash_nd, mut_nd;
	input_pt_inf(pt_inf);
	std::unordered_map<int, double> mut_p;
	std::vector<bool> prn;
	add_node(nd_muts, hash_nd, mut_nd, nd_par, nd_chd, nd_cnt, empty_vector);
	generate_nodes(pt_inf, nd_pt, nd_muts, hash_nd, mut_nd, nd_par, nd_chd, nd_cnt);
	for(std::unordered_map<int, std::vector<int>>::iterator i = mut_nd.begin(); i != mut_nd.end(); ++i) {
		int sum = 0; // sum = i.second[j].size(); weight by nodes
		for(int j = 0; j < i->second.size(); ++j)
			sum += (int) nd_pt[i->second[j]].size(); // weight by patient cases
		mut_p[i->first] = (double) sum / pt_inf.size();
	}
	prune_nodes(nd_muts, hash_nd, mut_nd, nd_par, nd_chd, nd_cnt, prn);
	back_track(nd_muts, hash_nd, mut_nd, nd_par, nd_chd, nd_cnt, prn, mut_p);
	output_parent(pt_inf, nd_pt, nd_muts, nd_par, nd_cnt);
	// output newick tree
	std::fstream newick_file;
	newick_file.open("newick_file.txt", std::fstream::out);
	dfs_newick(newick_file, pt_inf, nd_pt, nd_muts, nd_chd, nd_cnt, mut_p, 0);
	newick_file << ';';
	newick_file.close();
	// output text tree (directory like)
	std::fstream txt_file;
	txt_file.open("hCoV-2019 phylogenic tree.txt", std::fstream::out);
	dfs_txt(txt_file, pt_inf, nd_pt, nd_muts, nd_par, nd_chd, nd_cnt, mut_p, 0, std::string());
	txt_file.close();
	return 0;
}