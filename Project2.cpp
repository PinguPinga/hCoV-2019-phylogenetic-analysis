#define _CRT_SECURE_NO_WARNINGS
#include <utility>
#include <string>
#include <cmath>
#include <vector>
#include <iomanip>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <sstream>
#define H 196613

int thres_nodes, thres_leaves;
bool uniq_mut, show_all, show_cnt, show_warn;
bool output_sort_by_time, kind_taiwan;

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
	double p = 0; // p = 1;
	for(int i = 0; i < muts.size(); ++i)
		p += mut_p[muts[i]]; // p *= mut_p[muts[i]];
	return p;
}

std::string muts_str(std::vector<int> &muts) {
	if(muts.empty())
		return "NIL";
	std::string str;
	for(int i = 0; i < muts.size(); ++i) {
		str += num_mut(muts[i]);
		if(i != muts.size() - 1)
			str += ' ';
	}
	return str;
}

bool input_pt_inf(std::vector<inf> &pt_inf) {
	std::fstream input_file;
	input_file.open("Project2_input_file.txt", std::fstream::in);
	if(!input_file.is_open()) {
		std::cout << "input data reading error" << std::endl;
		return 0;
	}
	int n = 0, muts;
	// std::set<std::string> blacklist = {""};
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
		for(int i = 0; i < pt_inf[n].mut.size(); ++i) {
			input_file >> pt_inf[n].mut[i];
			if(pt_inf[n].mut[i] == 3863) { // reject all C241T
				--i;
				pt_inf[n].mut.pop_back();
			}
		}
		if(pt_inf[n].area == "belgium" || pt_inf[n].area == "netherlands") {
			--n;
			pt_inf.pop_back();
		}
		++n;
	}
	input_file.close();
	return 1;
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
	if(nd >= prn.size())
		prn.resize(nd + 1, 0);
	if(prn[nd] == 1)
		return;
	prn[nd] = 1;
	nd_cnt[nd] = 0;
}


void prune_node(std::vector<std::vector<int>> &nd_muts, 
	std::unordered_map<int, std::vector<int>> &hash_nd, 
	std::unordered_map<int, std::vector<int>> &mut_nd, 
	std::vector<int> &nd_par, 
	std::vector<std::unordered_set<int>> &nd_chd, 
	std::vector<int> &nd_cnt, 
	std::vector<bool> &prn, int nd) {
	if(prn[nd] == 1)
		return;
	prn[nd] = 1;
	std::vector<int> tmp_muts = nd_muts[nd];
	for(auto j : nd_muts[nd])
		if(mut_nd[j].size() == 1)
			tmp_muts.erase(std::find(tmp_muts.begin(), tmp_muts.end(), j));
	add_node(nd_muts, hash_nd, mut_nd, nd_par, nd_chd, nd_cnt, tmp_muts);
	int tmp_nd = find_nd(nd_muts, hash_nd, tmp_muts);
	nd_chd[tmp_nd].insert(nd);
	if(nd_cnt[tmp_nd] > 0) {
		nd_cnt[tmp_nd] += nd_cnt[nd];
		nd_par[nd] = tmp_nd;
	}
	// if(tmp_nd < prn.size())
	// 	prn[tmp_nd] = 0;
}

void prune_nodes(std::vector<inf> &pt_inf, 
	std::vector<std::vector<int>> &nd_pt, 
	std::vector<std::vector<int>> &nd_muts, 
	std::unordered_map<int, std::vector<int>> &hash_nd, 
	std::unordered_map<int, std::vector<int>> &mut_nd, 
	std::vector<int> &nd_par, 
	std::vector<std::unordered_set<int>> &nd_chd, 
	std::vector<int> &nd_cnt, 
	std::vector<bool> &prn) {
	prn.resize(nd_muts.size(), 0);
	if(uniq_mut)
		for(std::unordered_map<int, std::vector<int>>::iterator i = mut_nd.begin(); i != mut_nd.end(); ++i)
			if(i->second.size() == 1) { // del node with this unique mutation
				if(kind_taiwan) {
					bool pass = 0;
					for(auto j : nd_pt[i->second[0]])
						if(pt_inf[j].area == "taiwan") {
							pass = 1;
							continue;
						}
					if(pass)
						continue;
				}
				if(nd_cnt[i->second[0]] < thres_leaves)
					del_node(nd_muts, hash_nd, mut_nd, nd_par, nd_chd, nd_cnt, prn, i->second[0]);
				else
					prune_node(nd_muts, hash_nd, mut_nd, nd_par, nd_chd, nd_cnt, prn, i->second[0]);
			}
	if(thres_nodes > 1 || thres_leaves > 1)
		for(int i = 0; i < nd_cnt.size(); ++i) {
			if(kind_taiwan) {
				bool pass = 0;
				if(i < nd_pt.size())
					for(auto j : nd_pt[i])
						if(pt_inf[j].area == "taiwan") {
							pass = 1;
							continue;
						}
				if(pass)
					continue;
				for(auto j : nd_chd[i])
					for(auto k : nd_pt[j])
						if(pt_inf[k].area == "taiwan") {
							pass = 1;
							continue;
						}
				if(pass)
					continue;
			}
			int sum = nd_cnt[i];
			for(auto j : nd_chd[i])
				sum += nd_cnt[j];
			if(sum < thres_nodes && nd_cnt[i] < thres_leaves)
				del_node(nd_muts, hash_nd, mut_nd, nd_par, nd_chd, nd_cnt, prn, i);
		}
}

void back_track(std::vector<inf> &pt_inf, 
	std::vector<std::vector<int>> &nd_pt, 
	std::vector<std::vector<int>> &nd_muts, 
	std::unordered_map<int, std::vector<int>> &hash_nd, 
	std::unordered_map<int, std::vector<int>> &mut_nd, 
	std::vector<int> &nd_par, 
	std::vector<std::unordered_set<int>> &nd_chd, 
	std::vector<int> &nd_cnt, 
	std::vector<bool> &prn, 
	std::unordered_map<int, double> &mut_p) {
	std::vector<bool> inq(nd_muts.size());
	std::vector<std::vector<int>> nd_q;
	for(int i = 0; i < nd_muts.size(); ++i) { // node 0 is root {no mutation}
		if(nd_par[i] == -1 && nd_muts[i].size() != 0 && (i >= prn.size() || prn[i] == 0)) {
			if(nd_q.size() < nd_muts[i].size() + 1)
				nd_q.resize(nd_muts[i].size() + 1);
			nd_q[nd_muts[i].size()].emplace_back(i);
			inq[i] = 1;
		}
	}
	for(int i = (int) nd_q.size() - 1; i > 0; --i) { // layer of queue = mutations
		sort(nd_q[i].begin(), nd_q[i].end(), [&](int &a, int &b) {
			if(a < nd_pt.size() && b >= nd_pt.size())
				return true;
			else if(a >= nd_pt.size() && b < nd_pt.size())
				return false;
			else if(a < nd_pt.size() && b < nd_pt.size()) {
				if(pt_inf[nd_pt[a][0]].y != pt_inf[nd_pt[b][0]].y)
					return (pt_inf[nd_pt[a][0]].y < pt_inf[nd_pt[b][0]].y);
				if(pt_inf[nd_pt[a][0]].m != pt_inf[nd_pt[b][0]].m)
					return (pt_inf[nd_pt[a][0]].m < pt_inf[nd_pt[b][0]].m);
				if(pt_inf[nd_pt[a][0]].d != pt_inf[nd_pt[b][0]].d)
					return (pt_inf[nd_pt[a][0]].d < pt_inf[nd_pt[b][0]].d);
			}
			double a_p = muts_p(mut_p, nd_muts[a]), b_p = muts_p(mut_p, nd_muts[b]);
			if(a_p != b_p)
				return a_p > b_p;
			return a < b;
		} );
		for(int j = 0; j < nd_q[i].size(); ++j) {
			int nd = nd_q[i][j];
			for(std::unordered_set<int>::iterator k = nd_chd[nd].begin(); k != nd_chd[nd].end();) // refresh child list
				if(nd_par[*k] != -1 && nd_par[*k] != nd)
					nd_chd[nd].erase(k++);
				else
					++k;
			if(nd_chd[nd].empty() && nd_cnt[nd] == 0) // no need to trace back if no child and non-existent
				continue;
			if(nd_chd[nd].size() > 1 && nd_cnt[nd] == 0) { // create common ancestor of >=2 nodes
				if(nd < nd_pt.size())
					nd_cnt[nd] += nd_pt[nd].size();
				for(auto k : nd_chd[nd]) {
					nd_par[k] = nd;
					nd_cnt[nd] += nd_cnt[k];
				}
			}
			int lnk_nd = nd_cnt[nd] > 0 ? nd : *nd_chd[nd].begin();
			if(nd_par[lnk_nd] != -1)
				continue;
			std::vector<int> exst_lst, nexst_lst;
			int tmp_nd;
			std::vector<int> tmp_muts(nd_muts[nd].size() - 1);
			for(int k = 0; k + 1 < nd_muts[nd].size(); ++k) // copy nd_muts[1 .. last]
				tmp_muts[k] = nd_muts[nd][k + 1];
			add_node(nd_muts, hash_nd, mut_nd, nd_par, nd_chd, nd_cnt, tmp_muts);
			tmp_nd = find_nd(nd_muts, hash_nd, tmp_muts);
			if(tmp_nd >= nd_cnt.size() || (tmp_nd < nd_cnt.size() && nd_cnt[tmp_nd] == 0))
				nexst_lst.emplace_back(tmp_nd);
			else
				exst_lst.emplace_back(tmp_nd);
			for(int k = 0; k + 1 < nd_muts[nd].size(); ++k) {
				tmp_muts[k] = nd_muts[nd][k];
				add_node(nd_muts, hash_nd, mut_nd, nd_par, nd_chd, nd_cnt, tmp_muts);
				tmp_nd = find_nd(nd_muts, hash_nd, tmp_muts);
				if(tmp_nd >= nd_cnt.size() || (tmp_nd < nd_cnt.size() && nd_cnt[tmp_nd] == 0))
					nexst_lst.emplace_back(tmp_nd);
				else
					exst_lst.emplace_back(tmp_nd);
			}
			if(exst_lst.size() > 0) {
				tmp_nd = *min_element(exst_lst.begin(), exst_lst.end(), [&](int &a, int &b) {
				if(a < nd_pt.size() && b >= nd_pt.size())
					return true;
				else if(a >= nd_pt.size() && b < nd_pt.size())
					return false;
				else if(a < nd_pt.size() && b < nd_pt.size()) {
					if(pt_inf[nd_pt[a][0]].y != pt_inf[nd_pt[b][0]].y)
						return (pt_inf[nd_pt[a][0]].y < pt_inf[nd_pt[b][0]].y);
					if(pt_inf[nd_pt[a][0]].m != pt_inf[nd_pt[b][0]].m)
						return (pt_inf[nd_pt[a][0]].m < pt_inf[nd_pt[b][0]].m);
					if(pt_inf[nd_pt[a][0]].d != pt_inf[nd_pt[b][0]].d)
						return (pt_inf[nd_pt[a][0]].d < pt_inf[nd_pt[b][0]].d);
				}
				double a_p = muts_p(mut_p, nd_muts[a]), b_p = muts_p(mut_p, nd_muts[b]);
				if(a_p != b_p)
					return a_p > b_p;
				return a < b;
				} );
				nd_chd[tmp_nd].insert(lnk_nd);
				nd_par[lnk_nd] = tmp_nd;
				nd_cnt[tmp_nd] += nd_cnt[lnk_nd];
				if(tmp_nd >= inq.size())
					inq.resize(tmp_nd + 1, 0);
				if(inq[tmp_nd] == 0) {
					inq[tmp_nd] = 1;
					nd_q[tmp_muts.size()].emplace_back(tmp_nd);
				}
				if(exst_lst.size() > 1) {
					printf("node #%d has difficulty in determining parent. \n", nd);
					for(int k = 0; k < exst_lst.size(); ++k) {
						printf("%d", exst_lst[k]);
						if(k != exst_lst.size() - 1)
							printf("	");
						else
							printf("\n");
					}
				}
			}
			else
				for(int k = 0; k < nexst_lst.size(); ++k) {
					nd_chd[nexst_lst[k]].insert(lnk_nd);
					if(nexst_lst[k] >= inq.size())
						inq.resize(nexst_lst[k] + 1, 0);
					if(inq[nexst_lst[k]] == 0) {
						inq[nexst_lst[k]] = 1;
						nd_q[tmp_muts.size()].emplace_back(nexst_lst[k]);
					}
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
	std::fstream out_file;
	out_file.open("hCoV-2019 phylogenic tree parent data.txt", std::fstream::out);
	for(int i = 0; i < nd_muts.size(); ++i) {
		if(nd_cnt[i] == 0)
			continue;
		if(i < nd_pt.size() && nd_par[i] != -1 && nd_muts[i].size() - nd_muts[nd_par[i]].size() > 2) // nd_pt > 1 ? or subtree > 1?
			out_file << "WARNING: multiple (>=3) additional mutations" << std::endl;
		if(i >= nd_pt.size())
			out_file << "generated ";
		out_file << "node #" << i << std::endl;
		if(i < nd_pt.size())
			for(int j = 0; j < nd_pt[i].size(); ++j)
				out_file << pt_inf[nd_pt[i][j]].acc << ' ' << pt_inf[nd_pt[i][j]].area << ' '<< pt_inf[nd_pt[i][j]].y << '-'<< std::setfill('0') << std::setw(2) << pt_inf[nd_pt[i][j]].m << '-'<< std::setfill('0') << std::setw(2) << pt_inf[nd_pt[i][j]].d << std::endl;
		out_file << muts_str(nd_muts[i]) << std::endl;
		if(nd_par[i] >= nd_pt.size())
			out_file << "generated ";
		out_file << "parent #" << nd_par[i] << std::endl;
		if(nd_par[i] == -1)
			out_file << "NIL" << std::endl;
		else
			out_file << muts_str(nd_muts[nd_par[i]]) << std::endl;
		if(i != nd_muts.size() - 1)
			out_file << std::endl;
	}
	out_file.close();
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

void output_newick_tree(std::vector<inf> &pt_inf, 
	std::vector<std::vector<int>> &nd_pt, 
	std::vector<std::vector<int>> &nd_muts, 
	std::vector<std::unordered_set<int>> &nd_chd, 
	std::vector<int> &nd_cnt, 
	std::unordered_map<int, double> &mut_p) {
	std::fstream out_file;
	out_file.open("hCoV-2019 newick tree.txt", std::fstream::out);
	dfs_newick(out_file, pt_inf, nd_pt, nd_muts, nd_chd, nd_cnt, mut_p, 0);
	out_file << ';';
	out_file.close();
}

void dfs_phy(std::fstream &out_file, 
	std::vector<inf> &pt_inf, 
	std::vector<std::vector<int>> &nd_pt, 
	std::vector<std::vector<int>> &nd_muts, 
	std::vector<int> &nd_par, 
	std::vector<std::unordered_set<int>> &nd_chd, 
	std::vector<int> &nd_cnt, 
	std::vector<std::vector<int>> &nd_add, 
	std::unordered_map<int, double> &mut_p, 
	int nd, 
	std::string prefix, int prefix_len) {
	std::stringstream ss;
	ss << "node " << "#" << nd;
	if(nd != 0) { // output additional mutations
		ss << " (+" << muts_str(nd_add[nd]) << ")";
	}
	if(nd < nd_pt.size()) { // output patient data
		while(ss.str().length() + prefix_len < 40 + 20)
			ss << ' ';
		int t = 0;
		ss << "  " << pt_inf[nd_pt[nd][0]].acc << "  " << pt_inf[nd_pt[nd][0]].area;
		while(ss.str().length() + prefix_len < 71 + 20)
			ss << ' ';
		ss << ' ' << pt_inf[nd_pt[nd][0]].y << '-'<< std::setfill('0') << std::setw(2) << pt_inf[nd_pt[nd][0]].m << '-' << std::setfill('0') << std::setw(2) << pt_inf[nd_pt[nd][0]].d;
		if(show_cnt && (nd_pt[nd].size() > 1 || nd_cnt[nd] > nd_pt[nd].size() + nd_chd[nd].size())) {
			ss << " (";
			if(nd_pt[nd].size() == 1)
				ss << "1 entry"; 
			else
				ss << nd_pt[nd].size() << " entries"; 
			if(nd_cnt[nd] > nd_pt[nd].size() + nd_chd[nd].size())
				ss << ", " << nd_cnt[nd] << " entries under clade";
			ss << ")";
		}
		if(show_warn && (nd_par[nd] != -1 && nd_muts[nd].size() - nd_muts[nd_par[nd]].size() > 2)) { // output WARNING
			while(ss.str().length() + prefix_len < 200)
				ss << ' ';
			ss << "	WARNING: multiple (>=3) additional mutations";
		}
		out_file << ss.str() << std::endl;
		if(show_all) {
			for(int i = 1; i < nd_pt[nd].size(); ++i) {
				std::stringstream tmp_ss;
				out_file << prefix;
				int add_len = 0;
				if(!nd_chd[nd].empty()) {
					out_file << "¢x";
					add_len = 1;
				}
				while(tmp_ss.str().length() + prefix_len + add_len < 40 + 20)
					tmp_ss << ' ';
				int t = 0;
				tmp_ss << "  " << pt_inf[nd_pt[nd][i]].acc << "  " << pt_inf[nd_pt[nd][i]].area;
				while(tmp_ss.str().length() + prefix_len + add_len < 71 + 20)
					tmp_ss << ' ';
				tmp_ss << ' ' << pt_inf[nd_pt[nd][i]].y << '-'<< std::setfill('0') << std::setw(2) << pt_inf[nd_pt[nd][i]].m << '-' << std::setfill('0') << std::setw(2) << pt_inf[nd_pt[nd][i]].d;
				out_file << tmp_ss.str() << std::endl;
				tmp_ss.str("");
			}
		}
	}
	else {
		/* if(nd_par[nd] != -1 && nd_muts[nd].size() - nd_muts[nd_par[nd]].size() > 2) { // output WARNING
		while(ss.str().length() + prefix_len < 200)
			ss << ' ';
		ss << "	WARNING: multiprefix_lene (>=3) additional mutations";
		} */
		out_file << ss.str() << std::endl;
	}
	ss.str("");
	if(nd_chd[nd].empty())
		return;
	std::vector<int> lst;
	for(std::unordered_set<int>::iterator i = nd_chd[nd].begin(); i != nd_chd[nd].end(); ++i)
		lst.emplace_back(*i);
	sort(lst.begin(), lst.end(), [&] (int a, int b) {
		if(output_sort_by_time){ // sort by first case time
			while(a >= nd_pt.size() || nd_pt[a].empty())
				a = *nd_chd[a].begin();
			while(b >= nd_pt.size() || nd_pt[b].empty())
				b = *nd_chd[b].begin();
			if(pt_inf[nd_pt[a][0]].y != pt_inf[nd_pt[b][0]].y)
				return pt_inf[nd_pt[a][0]].y < pt_inf[nd_pt[b][0]].y;
			else if(pt_inf[nd_pt[a][0]].m != pt_inf[nd_pt[b][0]].m)
				return pt_inf[nd_pt[a][0]].m < pt_inf[nd_pt[b][0]].m;
			else if(pt_inf[nd_pt[a][0]].d != pt_inf[nd_pt[b][0]].d)
				return pt_inf[nd_pt[a][0]].d < pt_inf[nd_pt[b][0]].d;
		}
		if(nd_cnt[a] != nd_cnt[b])
			return nd_cnt[a] > nd_cnt[b];
		double a_p = muts_p(mut_p, nd_muts[a]), b_p = muts_p(mut_p, nd_muts[b]);
		if(a_p != b_p)
			return a_p < b_p;
		return a < b;
	} );
	std::string tmp = prefix;
	for(int i = 0; i < lst.size(); ++i) {
		out_file << prefix;
		tmp = prefix;
		int tt = prefix_len;
		if(i != lst.size() - 1) {
			out_file << "¢u";
			tmp += "¢x";
			++tt;
		}
		else if(i == lst.size() - 1) {
			out_file << "¢|";
			tmp += " ";
			++tt;
		}
		for(int j = 0; j < nd_muts[lst[i]].size() - (prefix.length() / 2) - 1; ++j) {
			out_file << "¢w";
			tmp += " ";
			++tt;
		}
		dfs_phy(out_file, pt_inf, nd_pt, nd_muts, nd_par, nd_chd, nd_cnt, nd_add, mut_p, lst[i], tmp, tt);
	}
}

void output_phy_tree(std::vector<inf> &pt_inf, 
	std::vector<std::vector<int>> &nd_pt, 
	std::vector<std::vector<int>> &nd_muts, 
	std::vector<int> &nd_par, 
	std::vector<std::unordered_set<int>> &nd_chd, 
	std::vector<int> &nd_cnt, 
	std::vector<std::vector<int>> &nd_add, 
	std::unordered_map<int, double> &mut_p) {
	// std::fstream out_file;
	// out_file.open("hCoV-2019 phylogenic tree.txt", std::fstream::out);
	// dfs_phy(out_file, pt_inf, nd_pt, nd_muts, nd_par, nd_chd, nd_cnt, nd_add, mut_p, 0, std::string(), 0);
	// out_file.close();
	std::fstream simp_file, comp_file;
	simp_file.open("hCoV-2019 phylogenic tree (simple).txt", std::fstream::out);
	show_all = 0;
	dfs_phy(simp_file, pt_inf, nd_pt, nd_muts, nd_par, nd_chd, nd_cnt, nd_add, mut_p, 0, std::string(), 0);
	simp_file.close();
	
	comp_file.open("hCoV-2019 phylogenic tree (complete).txt", std::fstream::out);
	show_all = 1;
	dfs_phy(comp_file, pt_inf, nd_pt, nd_muts, nd_par, nd_chd, nd_cnt, nd_add, mut_p, 0, std::string(), 0);
	comp_file.close();
}

void mut_p_init(std::unordered_map<int, std::vector<int>> &mut_nd, 
	std::vector<std::vector<int>> &nd_pt, 
	int n, 
	std::unordered_map<int, double> &mut_p) {
	for(std::unordered_map<int, std::vector<int>>::iterator i = mut_nd.begin(); i != mut_nd.end(); ++i) {
		int sum = 0; // sum = i.second[j].size(); weight by nodes
		for(int j = 0; j < i->second.size(); ++j)
			sum += (int) nd_pt[i->second[j]].size(); // weight by patient cases
		mut_p[i->first] = (double) n / sum; // (double) sum / n;
	}
}

void sort_node(std::vector<inf> &pt_inf, std::vector<std::vector<int>> &nd_pt) {
	for(int i = 0; i < nd_pt.size(); ++i) {
		std::unordered_map<std::string, int> area_cnt;
		for(int j = 0; j < nd_pt[i].size(); ++j)
			++area_cnt[pt_inf[nd_pt[i][j]].area];
		std::sort(nd_pt[i].begin(), nd_pt[i].end(), [&] (int &a, int &b) {
			if(pt_inf[a].y != pt_inf[b].y)
				return pt_inf[a].y < pt_inf[b].y;
			if(pt_inf[a].m != pt_inf[b].m)
				return pt_inf[a].m < pt_inf[b].m;
			if(pt_inf[a].d != pt_inf[b].d)
				return pt_inf[a].d < pt_inf[b].d;
			if(area_cnt[pt_inf[a].area] != area_cnt[pt_inf[b].area])
				return area_cnt[pt_inf[a].area] > area_cnt[pt_inf[b].area];
			int scmp = pt_inf[a].area.compare(pt_inf[b].area);
			if(scmp != 0)
				return scmp < 0;
			return a < b;
		} );
	}
}

void muts_diff(std::vector<std::vector<int>> &nd_muts, int a, int b, std::vector<int> &tmp_muts) { // a - b
	if(a < 0 || a >= nd_muts.size() || b < 0 || b >= nd_muts.size())
		return;
	int i = 0, j = 0;
	while(i < nd_muts[a].size() && j < nd_muts[b].size()) {
		if(nd_muts[a][i] < nd_muts[b][j])
			tmp_muts.emplace_back(nd_muts[a][i++]);
		else if(nd_muts[a][i] == nd_muts[b][j])
			++i, ++j;
		else if(nd_muts[a][i] > nd_muts[b][j])
			++j;
	}
	while(i < nd_muts[a].size())
		tmp_muts.emplace_back(nd_muts[a][i++]);
}

void add_init(std::vector<std::vector<int>> &nd_muts, std::vector<int> &nd_par, std::vector<std::vector<int>> &nd_add, std::unordered_map<int, std::vector<int>> &add_nd) {
	nd_add.resize(nd_par.size());
	for(int i = 0; i < nd_add.size(); ++i) {
		if(nd_par[i] != -1)
			muts_diff(nd_muts, i, nd_par[i], nd_add[i]);
		for(int j = 0; j < nd_add[i].size(); ++j)
			add_nd[nd_add[i][j]].emplace_back(i);
	}
}

void output_new_input(std::vector<inf> &pt_inf, 
	std::vector<std::vector<int>> &nd_pt, 
	std::vector<std::vector<int>> &nd_muts, 
	std::unordered_map<int, std::vector<int>> &add_nd) {
	std::fstream out_file;
	out_file.open("new input_file.txt", std::fstream::out);
	for(int i = 0; i < nd_pt.size(); ++i) {
		for(int j = 0; j < nd_pt[i].size(); ++j) {
			out_file << pt_inf[nd_pt[i][j]].acc << '	';
			out_file << pt_inf[nd_pt[i][j]].area << '	';
			out_file << pt_inf[nd_pt[i][j]].y << '	';
			out_file << pt_inf[nd_pt[i][j]].m << '	';
			out_file << pt_inf[nd_pt[i][j]].d << std::endl;
			std::stringstream ss;
			int n = 0;
			for(int k = 0; k < nd_muts[i].size(); ++k) {
				// if(add_nd.find(nd_muts[i][k]) == add_nd.end())
				//	continue;
				if(add_nd[nd_muts[i][k]].size() >= 2) {
					std::cout << num_mut(nd_muts[i][k]) << std::endl;
				}
				if(add_nd[nd_muts[i][k]].size() > 3) {
					continue;
				}
				ss << '	' << nd_muts[i][k];
				++n;
			}
			out_file << n << ss.str();
			out_file << std::endl << std::endl;
		}
	}
	system("pause");
	out_file.close();
}

void dfs_tree(
	std::vector<std::unordered_set<int>> &nd_chd, 
	int nd, std::vector<int> &nd_lst, std::vector<std::pair<int, int>> &edge_lst) {
	nd_lst.emplace_back(nd);
	for(auto i : nd_chd[nd]) {
		edge_lst.emplace_back(std::pair<int, int>(nd, i));
		dfs_tree(nd_chd, i, nd_lst, edge_lst);
	}
}

void draw_tree(std::vector<std::unordered_set<int>> &nd_chd) {
	std::fstream out_file;
	std::vector<int> nd_lst;
	std::vector<std::pair<int, int>> edge_lst;
	dfs_tree(nd_chd, 0, nd_lst, edge_lst);
	out_file.open("draw tree.txt", std::fstream::out);
	std::unordered_map<int, int> m;
	for(int i = 0; i < nd_lst.size(); ++i) {
		out_file << i << std::endl;
		m[nd_lst[i]] = i;
	}
	for(int i = 0; i < edge_lst.size(); ++i)
		out_file << m[edge_lst[i].first] << ' ' << m[edge_lst[i].second] << std::endl;
	out_file.close();
}

bool input_setting() {
	std::fstream setting_file;
	setting_file.open("Project2_setting_file.txt", std::fstream::in);
	if(!setting_file.is_open()) {
		std::cout << "no setting file" << std::endl;
		return 0;
	}
	std::string buf;
	std::stringstream ss;
	getline(setting_file, buf);
	ss.str(buf);
	ss >> thres_nodes; // filter nodes with less than < N cases (int)
	getline(setting_file, buf);
	ss.str(buf);
	ss >> thres_leaves; // filter nodes with less than < N cases (int)
	getline(setting_file, buf);
	ss.str(buf);
	ss >> uniq_mut; // prune nodes with unique mutations
	getline(setting_file, buf);
	ss.str(buf);
	ss >> show_all; // show all cases in phylogenetic tree
	getline(setting_file, buf);
	ss.str(buf);
	ss >> show_cnt; // show entry counts
	getline(setting_file, buf);
	ss.str(buf);
	ss >> show_warn; // show warnings
	getline(setting_file, buf);
	ss.str(buf);
	ss >> output_sort_by_time; // output sort by time, otherwise sort by clade possibility (calculated penalty)
	getline(setting_file, buf);
	ss.str(buf);
	ss >> kind_taiwan; // output sort by time, otherwise sort by clade possibility (calculated penalty)
	setting_file.close();
	return 1;
}

int main()
{
	std::vector<inf> pt_inf; // patient information
	std::vector<int> nd_par, nd_cnt, empty_vector; // node=>parent, node=>node counts in branch
	std::vector<std::vector<int>> nd_muts, nd_pt, nd_add; // node=>mutations, nd=>patients, nd_add
	std::vector<std::unordered_set<int>> nd_chd; // node=>children
	std::unordered_map<int, std::vector<int>> hash_nd, mut_nd, add_nd; // hash mutations=>node, mutation=>nodes which have this mutation
	std::unordered_map<int, double> mut_p; // mutation probability
	std::vector<bool> prn; // pruned?
	if(input_setting() == 0) { // input
		system("pause");
		return 0;
	}
	if(input_pt_inf(pt_inf) == 0) {
		system("pause");
		return 0;
	}
	add_node(nd_muts, hash_nd, mut_nd, nd_par, nd_chd, nd_cnt, empty_vector);
	generate_nodes(pt_inf, nd_pt, nd_muts, hash_nd, mut_nd, nd_par, nd_chd, nd_cnt);
	sort_node(pt_inf, nd_pt);
	mut_p_init(mut_nd, nd_pt, (int) pt_inf.size(), mut_p);
	prune_nodes(pt_inf, nd_pt, nd_muts, hash_nd, mut_nd, nd_par, nd_chd, nd_cnt, prn);
	back_track(pt_inf, nd_pt, nd_muts, hash_nd, mut_nd, nd_par, nd_chd, nd_cnt, prn, mut_p);
	add_init(nd_muts, nd_par, nd_add, add_nd);
	output_parent(pt_inf, nd_pt, nd_muts, nd_par, nd_cnt);
	//output_newick_tree(pt_inf, nd_pt, nd_muts, nd_chd, nd_cnt, mut_p); // output newick tree
	output_phy_tree(pt_inf, nd_pt, nd_muts, nd_par, nd_chd, nd_cnt, nd_add, mut_p);	// output phylogenetic text tree (directory like)
	//draw_tree(nd_chd);
	//output_new_input(pt_inf, nd_pt, nd_muts, add_nd);
	return 0;
}