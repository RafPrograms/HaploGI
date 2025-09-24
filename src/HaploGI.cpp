 /**
*    The HaploGI (Haplotyping Given Inheritance) program is designed to
*    perform pedigree-based haplotyping of WGS (Whole Genome Sequence) data
*    and determination of haplotype sharing among cases with evidence for linkage
*    as described in Nafikov et al. (2025).
*    Copyright (C) 2025  Rafael A. Nafikov  email <nrafscience@gmail.com>
*
*    This program is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/


#include <time.h>
#include <sstream>
#include <fstream>
#include <map>
#include <list>
#include <memory>
#include <cstdlib>
#include <string.h>
#include <ctime>
#include <regex>
#include <iterator>
#include <iomanip>
#include <algorithm>
#include <tuple>
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <string>
#include <math.h>
#include <unistd.h>
#include <thread>
#include <streambuf>
#include <set>
#include <cstring> 
#include <random> // for std::mt19937
#include "../../../../usr/include/linux/limits.h"

using namespace std;

/// this section finds difference, union, and intersection between two arrays (vectors)
template <class T1>
class Comp_arrays_venn
{

public:
    Comp_arrays_venn(std::vector<T1> *, std::vector<T1> *);

    std::vector<T1> getDifference();

    std::vector<T1> getUnion();

    std::vector<T1> getIntersection();

private:
    std::vector<T1> *array1 = nullptr;
    std::vector<T1> *array2 = nullptr;
    std::vector<T1> comb_array;
    std::vector<T1> comb_array_unique;
    std::vector<T1> intersection;
    std::vector<T1> difference;
};

template <class T1>
Comp_arrays_venn<T1>::Comp_arrays_venn(std::vector<T1> *ar1, std::vector<T1> *ar2) : array1{ar1}, array2{ar2}
{

    std::sort(array1->begin(), array1->end());
    std::sort(array2->begin(), array2->end());

    array1->erase(unique(array1->begin(), array1->end()), array1->end());
    array2->erase(unique(array2->begin(), array2->end()), array2->end());

    comb_array.insert(comb_array.begin(), array1->begin(), array1->end());
    comb_array.insert(comb_array.end(), array2->begin(), array2->end());
    std::sort(comb_array.begin(), comb_array.end());

    comb_array_unique.insert(comb_array_unique.begin(), comb_array.begin(), comb_array.end());
    comb_array_unique.erase(unique(comb_array_unique.begin(), comb_array_unique.end()), comb_array_unique.end());

    typename std::vector<T1>::iterator iter1;
    for (iter1 = comb_array_unique.begin(); iter1 != comb_array_unique.end(); ++iter1)
    {
        int ct1 = count(comb_array.begin(), comb_array.end(), *iter1);
        if (ct1 == 1)
        {

            difference.push_back(*iter1);
        }
        else
        {

            intersection.push_back(*iter1);
        }
    }
}

template <class T1>
std::vector<T1> Comp_arrays_venn<T1>::getDifference()
{

    std::vector<T1> result;

    result.insert(result.begin(), difference.begin(), difference.end());

    std::sort(result.begin(), result.end());

    return result;
}

template <class T1>
std::vector<T1> Comp_arrays_venn<T1>::getUnion()
{

    std::vector<T1> result;

    result.insert(result.begin(), comb_array_unique.begin(), comb_array_unique.end());

    std::sort(result.begin(), result.end());
    
    return result;
}

template <class T1>
std::vector<T1> Comp_arrays_venn<T1>::getIntersection()
{

    std::vector<T1> result;

    result.insert(result.begin(), intersection.begin(), intersection.end());

     std::sort(result.begin(), result.end());
    
    return result;
}


/// generates random numbers from the standard uniform distribution
double rand(std::mt19937& rng)
{
    std::uniform_real_distribution<double> u_distr(0.0, 1.0);
    double randNum = u_distr(rng);

    return randNum;
}


/// this section finds all possible combinations for a set of numbers
class Combinations
{

public:
    void getCombs(int, int);
    std::vector<std::vector<int>> *getAllCombs();
    Combinations(std::vector<int> *);

private:
    std::vector<int> *init_numbs = nullptr;
    std::vector<std::vector<int>> all_combs;
    std::vector<int> comb;
};

Combinations::Combinations(std::vector<int> *input_num) : init_numbs{input_num}
{

    std::sort(init_numbs->begin(), init_numbs->end());
}

void Combinations::getCombs(int offset_in, int k_in)
{
    int k = k_in;
    int offset = offset_in;

    if (k == 0)
    {
        all_combs.push_back(comb);
        return;
    }
    for (size_t i = offset; i <= init_numbs->size() - k; ++i)
    {
        comb.push_back((*init_numbs)[i]);
        getCombs(i + 1, k - 1);
        comb.pop_back();
    }
}

std::vector<std::vector<int>> *Combinations::getAllCombs()
{

    std::vector<std::vector<int>> *result = nullptr;

    int el_num = init_numbs->size();

    for (int k = 1; k <= el_num; k++)
    {

        this->getCombs(0, k);
    }

    result = &all_combs;
    return result;
}


double median(std::vector<int> &vectorIn)
{

    double result1;

    std::sort(vectorIn.begin(), vectorIn.end());

    if (vectorIn.size() % 2 == 0)
    {
        // even
        int temp1 = vectorIn.size() / 2;
        result1 = (vectorIn[temp1 - 1] + vectorIn[temp1]) / 2;
    }

    else

    {
        // odd
        int temp1 = ((vectorIn.size() - 1) / 2);
        result1 = vectorIn[temp1];
    }

    return result1;
}

typedef std::vector<std::string> GenoVec;
typedef std::map<int, std::map<int, GenoVec>> DenseGenoMap;

struct Pedigree
{ // sequential numbers from 0 to # of subjects in pedigree - 1
    std::string subject_id;
    std::string father_id;
    std::string mother_id;
    int subject_sex;
    int mat_mi;            // 0 = absent or actual number
    int pat_mi;            // 0 = absent or actual number
    int mat_fgl;           // fgl label for founder chr; 0 if not a founder chr
    int pat_fgl;           // fgl label for founder chr; 0 if not a founder chr
    int wgs_data = {};     // 0 = no, 1 = yes; absence or presence of WGS data
    int mat_3D_pos = {-1}; // indicates position in 3D matrix and only subjects with WGS data will have number other than -1
    int pat_3D_pos = {-1}; // indicates position in 3D matrix and only subjects with WGS data will have number other than -1
    int pheno = {-1};      // -1 = no phenotype data available, 0 = missing, 1 = unaffected, 2 = affected
    std::string subject_id_orig;// subject id as it appears in input files
    std::string father_id_orig;// father id as it appears in input files
    std::string mother_id_orig;// mother id as it appears in input files
};

struct Parameters
{
    std::string pedigree_file_path;
    std::string meiosis_file_path;
    std::string dense_marker_pos_file_path;
    std::string dense_marker_geno_txt_file_path;
    std::string framework_marker_pos_file_path;
    std::string haplotype_file_path;
    std::string candidate_grp_file_path;
    std::string output_dir_path;
    std::vector<double> linkage_region;
    double maxLODmarker;
    std::vector<std::string> candidate_group;
    int run_type; // 1 = haplotyping; 2 = haplosharing; 3 = full (both haplotyping and haplosharing)
    int num_iter;
    int window_size = 20;
    int chr;
    int seed = 1234;
    int iter4phasing = -1;// when = -1 then the program will determine itself which iteration to use for phasing
    double minShareLength = 0.25; // min length of FGL sharing haplotype as the percentage of ROI
    Parameters() {}
};

struct ListNode // structure for singular linked lists
{
    int value;
    struct ListNode *next;
};

class NumberList
{
private:
    ListNode *head;

public:
    NumberList() { head = nullptr; }

    ~NumberList();

    void appendNode(int);
    void insertNode(int);
    void deleteNode(int);
    void displayList() const;
    int countListNodes();
    
    bool checkListEmpty();

    ListNode *ptrNodeElement(int);
    

    ListNode *begin();
    ListNode *end();
};

bool NumberList::checkListEmpty()
{
    ListNode *nodePtr = nullptr;
    if (nodePtr == nullptr)
    {
        return true;
    }
    else
    {
        return false;
    }
}

int NumberList::countListNodes()
{
    ListNode *nodePtr;
    nodePtr = head;
    int length = 0;
    while (nodePtr)
    {
        length++;
        nodePtr = nodePtr->next;
    }
    return (length);
}

ListNode *NumberList::ptrNodeElement(int element) // element numbering starts from 0
{
    int node_n = countListNodes();
    ListNode *nodePtr;
    nodePtr = head;
    int length = 0;
    if (element >= node_n)
    {
        exit(0);
    }

    while (nodePtr && element != length)
    {
        length++;
        nodePtr = nodePtr->next;
    }

    return (nodePtr);
}

ListNode *NumberList::begin() { return head; } 

ListNode *NumberList::end() { return 0; } 

void NumberList::appendNode(int num)
{
    ListNode *newNode; 
    ListNode *nodePtr; 
   
    newNode = new ListNode;
    newNode->value = num;
    newNode->next = nullptr;
    
    if (!head)
        head = newNode;
    else 
    {
        
        nodePtr = head;
        
        while (nodePtr->next)
        nodePtr = nodePtr->next;
        
        nodePtr->next = newNode;
    }
}

void NumberList::displayList() const
{
    ListNode *nodePtr; 
    
    nodePtr = head;
    
    while (nodePtr)
    {
        
        std::cout << nodePtr->value << std::endl;
       
        nodePtr = nodePtr->next;
    }
}

void NumberList::insertNode(int num)
{
    ListNode *newNode;
    
    ListNode *nodePtr;
   
    ListNode *previousNode = nullptr; 
    
    newNode = new ListNode;
    newNode->value = num;
    
    if (!head)
    {
        head = newNode;
        newNode->next = nullptr;
    }
    else 
    {
        
        nodePtr = head;
        
        previousNode = nullptr;
        
        while (nodePtr != nullptr && nodePtr->value < num)
        {
            previousNode = nodePtr;
            nodePtr = nodePtr->next;
        }
        
        if (previousNode == nullptr)
        {
            head = newNode;
            newNode->next = nodePtr;
        }
        else 
        {
            previousNode->next = newNode;
            newNode->next = nodePtr;
        }
    }
}

void NumberList::deleteNode(int num)
{
    ListNode *nodePtr;
    
    ListNode *previousNode; 
    
    if (!head)
        return;
    
    if (head->value == num)
    {
        nodePtr = head->next;
        delete head;
        head = nodePtr;
    }
    else
    {
        
        nodePtr = head;
        
        while (nodePtr != nullptr && nodePtr->value != num)
        {
            previousNode = nodePtr;
            nodePtr = nodePtr->next;
        }
        
        
       
        if (nodePtr)
        {
            previousNode->next = nodePtr->next;
            delete nodePtr;
        }
    }
}

NumberList::~NumberList()
{
    ListNode *nodePtr;
    
    ListNode *nextNode; 
    
    nodePtr = head;
    
    while (nodePtr != nullptr)
    {
       
        nextNode = nodePtr->next;
       
        delete nodePtr;
        
        nodePtr = nextNode;
    }
}


class Pedigree_info
{

public:
    std::vector<Pedigree> cur_ped;

    Pedigree_info(std::vector<Pedigree> &ped1) : cur_ped{ped1} {}

    int getNumMeiosisLines();

    int getSubjectIndex(std::string); // gets element number for particular subject id from Pedigree struct

    int getNumSubjects();

    int getNumFounders();

    int isFounder(int);

    int getNumChildren(std::string);

    void genMI_FGLrelInfo();

    void printPedigree();

    void generateTempIDs();

    std::string convert2TempID(std::string);

    std::string convert2OrigID(std::string);

    int getMIlineNum(int, int);

    int getNumWGS();

    int getNumWGSlist(std::vector<std::string> &);

    std::vector<int> getCases();
};

void Pedigree_info::printPedigree()
{

    std::cout << "subject_id" << "  " << "father_id" << "  " << "mother_id"
    << "  " << "subject_sex" << "  " << "phenotype" << "  " << "mat_mi_line"
    << "  " << "pat_mi_line" << "  " << "founder_mat_fgl" << "  " << "founder_pat_fgl"
    << "  " << "WGS_data" << "  " << "mat_3D_pos" << "  " << "pat_3D_pos" 
    << "  " << "subject_id_orig" << "  " << "father_id_orig" << "  " << "mother_id_orig" << std::endl;

    for (size_t i = 0; i < cur_ped.size(); ++i)
    {
        std::cout << cur_ped[i].subject_id << "  " << cur_ped[i].father_id << "  " << cur_ped[i].mother_id << "  "
        << cur_ped[i].subject_sex << "  " << cur_ped[i].pheno << "  " << cur_ped[i].mat_mi << "  " 
        << cur_ped[i].pat_mi << "  " << cur_ped[i].mat_fgl << "  " << cur_ped[i].pat_fgl << "  " 
        << cur_ped[i].wgs_data << " " << cur_ped[i].mat_3D_pos << "  " << cur_ped[i].pat_3D_pos << "  " 
        << cur_ped[i].subject_id_orig << "  " << cur_ped[i].father_id_orig << "  " << cur_ped[i].mother_id_orig << std::endl;
    }
}

void Pedigree_info::generateTempIDs()
{

    for (size_t i = 0; i < cur_ped.size(); ++i)
    {    
        if (cur_ped[i].father_id_orig.compare("0") != 0) // different
        {
            for (size_t j = 0; j < cur_ped.size(); ++j)
            {
                if (cur_ped[i].father_id_orig.compare(cur_ped[j].subject_id_orig) == 0)
                {
                    cur_ped[i].father_id = cur_ped[j].subject_id; 
                }
            }
        }
        else
        {
            cur_ped[i].father_id = "0";
        }


        if (cur_ped[i].mother_id_orig.compare("0") != 0) // different
        {
            for (size_t j = 0; j < cur_ped.size(); ++j)
            {
                if (cur_ped[i].mother_id_orig.compare(cur_ped[j].subject_id_orig) == 0)
                {
                    cur_ped[i].mother_id = cur_ped[j].subject_id; 
                }
            }
        }
        else
        {
            cur_ped[i].mother_id = "0";
        }  
    }


}

std::string Pedigree_info::convert2TempID(std::string subject_id_orig)
{
    std::string subject_id_temp;

    for (size_t i = 0; i < cur_ped.size(); ++i)
    {
        if (subject_id_orig.compare (cur_ped[i].subject_id_orig) == 0)
        {
            subject_id_temp = cur_ped[i].subject_id;
            break;
        }

    }

    return subject_id_temp;
}

std::string Pedigree_info::convert2OrigID(std::string subject_id_temp)
{
    std::string subject_id_orig;

    for (size_t i = 0; i < cur_ped.size(); ++i)
    {
        if (subject_id_temp.compare (cur_ped[i].subject_id) == 0)
        {
            subject_id_orig = cur_ped[i].subject_id_orig;
            break;
        }

    }
    
    return subject_id_orig;
}

std::vector<int> Pedigree_info::getCases()
{

    std::vector<int> result1;

    for (size_t i = 0; i < cur_ped.size(); ++i)
    {

        if (cur_ped[i].pheno == 2) // leave only AD subjects
        {
            result1.push_back(i);
        }
    }

    return result1;
}

int Pedigree_info::getNumWGSlist(std::vector<std::string> &sub_list)
{

    int result1 = 0;

    for (size_t i = 0; i < sub_list.size(); ++i)
    {
        if (cur_ped[this->getSubjectIndex(sub_list[i])].wgs_data == 1)
        {
            result1++;
        }
    }

    return result1;
}

int Pedigree_info::getMIlineNum(int index, int chr_origin)
{
    if (chr_origin == 0)
    {
        return cur_ped[index].mat_mi;
    }
    else if (chr_origin == 1)
    {
        return cur_ped[index].pat_mi;
    }
    else
    {
        return -1;
    }
}

int Pedigree_info::getNumChildren(std::string id)
{

    int children_num = 0;

    for (size_t i = 0; i < cur_ped.size(); ++i)
    {
        if (cur_ped[i].father_id.compare(id) == 0 || cur_ped[i].mother_id.compare(id) == 0)
        {
            children_num++;
        }
    }

    return children_num;
}

void Pedigree_info::genMI_FGLrelInfo()
{

    int mi_num = 0;
    int fgl_num = 0;

    for (size_t i = 0; i < cur_ped.size(); ++i)
    {

        if (isFounder(i) == 1)
        {

            cur_ped[i].mat_mi = 0;
            fgl_num++;
            cur_ped[i].mat_fgl = fgl_num;

            cur_ped[i].pat_mi = 0;
            fgl_num++;
            cur_ped[i].pat_fgl = fgl_num;
        }

        else if (isFounder(i) == 0)
        {
            cur_ped[i].mat_fgl = 0;
            cur_ped[i].pat_fgl = 0;

            int mother_index = getSubjectIndex(cur_ped[i].mother_id);
            int father_index = getSubjectIndex(cur_ped[i].father_id);

            if (isFounder(mother_index) == 1 && getNumChildren(cur_ped[mother_index].subject_id) == 1)
            {
                mi_num++;
                cur_ped[i].mat_mi = 0;
                cur_ped[i].pat_mi = mi_num;
            }

            else if (isFounder(father_index) == 1 && getNumChildren(cur_ped[father_index].subject_id) == 1)
            {
                mi_num++;
                cur_ped[i].pat_mi = 0;
                cur_ped[i].mat_mi = mi_num;
            }

            else

            {
                mi_num++;
                cur_ped[i].mat_mi = mi_num;
                mi_num++;
                cur_ped[i].pat_mi = mi_num;
            }
        }
    }
}

int Pedigree_info::getNumFounders()
{

    int founders_num = 0;

    for (size_t i = 0; i < cur_ped.size(); ++i)
    {
        if (cur_ped[i].father_id.compare("0") == 0 && cur_ped[i].mother_id.compare("0") == 0)
        {
            founders_num++;
        }
    }

    return founders_num;
}

int Pedigree_info::getNumWGS()
{

    int wgs_num = 0;

    for (size_t i = 0; i < cur_ped.size(); ++i)
    {
        if (cur_ped[i].wgs_data == 1)
        {
            wgs_num++;
        }
    }

    return wgs_num;
}

int Pedigree_info::getSubjectIndex(std::string id)
{

    for (size_t i = 0; i < cur_ped.size(); ++i)
    {
        if (cur_ped[i].subject_id.compare(id) == 0)
        {
            return i;
        }
    }

    return -1;
}

int Pedigree_info::getNumMeiosisLines()
{

    int meiosis_lines_num = 0;
    for (size_t i = 0; i < cur_ped.size(); ++i)
    {
        if (isFounder(i) == 0)
        {
            int mother_index = getSubjectIndex(cur_ped[i].mother_id);
            int father_index = getSubjectIndex(cur_ped[i].father_id);

            if ((isFounder(mother_index) == 1 && getNumChildren(cur_ped[mother_index].subject_id) == 1) || (isFounder(father_index) == 1 && getNumChildren(cur_ped[father_index].subject_id) == 1))
            {
                meiosis_lines_num++;
            }

            else

            {
                meiosis_lines_num += 2;
            }
        }
    }
    return meiosis_lines_num;
}

int Pedigree_info::getNumSubjects()
{
    return cur_ped.size();
}

int Pedigree_info::isFounder(int index)
{

    if (cur_ped[index].father_id.compare("0") == 0 && cur_ped[index].mother_id.compare("0") == 0)
    {
        return (1); // yes
    }
    return (0); // no
}


class DenseMarkerPos
{

public:
    void setVal(int, double);

    void makeWindows();

    double getPos_cM(int);

    int getPos_bp_left(double);

    std::vector<int> getWindowBps(int);

    int getWindowWithBp(int);

    int getNumWindows();

    std::vector<int> getROIboundWindNum(std::vector<double>);

    DenseMarkerPos(int);

private:
    int window_s;
    std::map<int, double> dense_marker_pos;
    std::vector<std::vector<int>> windows;
};

DenseMarkerPos::DenseMarkerPos(int window_size) : window_s{window_size},
dense_marker_pos{}, windows{}
{
}

void DenseMarkerPos::makeWindows()
{

    std::map<int, double>::iterator iter1;

    int count1 = 0;

    int window_count = 0;

    std::vector<int> temp1;

    int windows_total = dense_marker_pos.size() / window_s;

    for (iter1 = dense_marker_pos.begin(); iter1 != dense_marker_pos.end(); iter1++)
    {
        
        if (window_count < windows_total)
        {
            count1++;
            if (count1 < window_s)
            {
                temp1.push_back(iter1->first);
            }

            else if (count1 == window_s)
            {
                temp1.push_back(iter1->first);
                windows.push_back(temp1);
                temp1.clear();
                count1 = 0;
                window_count++;
            }
        }
    }

    if (temp1.size() > 0)
    {
        windows.push_back(temp1);
    }
}

std::vector<int> DenseMarkerPos::getWindowBps(int window_num)
{

    std::vector<int> result;

    for (size_t i = 0; i < windows[window_num - 1].size(); i++)
    {
        result.push_back(windows[window_num - 1][i]);
    }

    return result;
}

int DenseMarkerPos::getNumWindows()
{

    return windows.size();
}

void DenseMarkerPos::setVal(int bp_pos_in, double cM_pos_in)
{

    dense_marker_pos.insert(std::make_pair(bp_pos_in, cM_pos_in));
}

double DenseMarkerPos::getPos_cM(int bp_pos_in)
{

    double cM_pos_out;

    cM_pos_out = dense_marker_pos.find(bp_pos_in)->second;

    return cM_pos_out;
}

/// returns first bp with target cM value or -1 if it is not found
int DenseMarkerPos::getPos_bp_left(double cM_pos_in)
{

    int bp_pos_out;

    std::vector<int> temp1;

    for (auto& it : dense_marker_pos) 
    {
        if (it.second == cM_pos_in)
        {
            temp1.push_back (it.first);
        }
    }

    sort(temp1.begin(), temp1.end()); 

    if (!temp1.empty())
    {
        bp_pos_out = temp1[0];
    }
    
    else
    
    {
       bp_pos_out = -1; 
    }                          
  
    return bp_pos_out;
}

/// gets window # (starting from 1) containing specified bp position or returns -1
int DenseMarkerPos::getWindowWithBp(int bp_pos_in)
{

    int result1 = -1;

    for(size_t i=0;i<windows.size();i++)
    {
		for(size_t j=0;j<windows[i].size();j++)
        {
            if (windows[i][j] == bp_pos_in)
            {
                result1 = i + 1;
                goto LB1;
            }
        }
    }

    LB1:;    

    return result1;
}


std::vector<int> DenseMarkerPos::getROIboundWindNum(std::vector<double> linkageRegion)
{

    std::vector<int> result1; // window numbers starts with 1
    int bpLeftBound, leftWnd, bpRightBound, rightWnd;

    bpLeftBound = this->getPos_bp_left(linkageRegion[0]);
    bpRightBound = this->getPos_bp_left(linkageRegion[1]);

    if ((bpLeftBound == -1) || (bpRightBound == -1))
    {
        result1 = {-1, -1};
        goto LB1;
    }

    leftWnd = this->getWindowWithBp(bpLeftBound);
    rightWnd = this->getWindowWithBp(bpRightBound);

    if ((leftWnd == -1) || (rightWnd == -1))
    {
        result1 = {-1, -1};
        goto LB1;
    }

    result1.push_back(leftWnd);
    result1.push_back(rightWnd);


    LB1:;
    return (result1); // returns -1s if it did not determine ROI boundaries in window #

}


class MI_matrix
{

public:
    void setVal(int, int, int, int, int);

    void printMeiosisLine(int, int);

    std::vector<int> getMI(int, int);

    void setRecombPos(int, int, int, std::vector<int>);

    void printRecombPos(int, int, int);

    std::vector<int> getRecombPos(int, int);

    MI_matrix() {}

    MI_matrix(int, int);

    ~MI_matrix();

private:
    int dimn1, dimn2;
    int ***matrixMI = nullptr;
    std::vector<std::vector<std::vector<int>>> recombPos;
};

MI_matrix::MI_matrix(int iter_num, int m_line_num) : dimn1{iter_num}, dimn2{m_line_num}
{

    matrixMI = new int **[dimn1];
    for (int i = 0; i < dimn1; i++)
    {
        matrixMI[i] = new int *[dimn2];
        for (int j = 0; j < dimn2; j++)
        {
            matrixMI[i][j] = new int[3];
            for (int k = 0; k < 3; k++)
            {
                matrixMI[i][j][k] = -1;
            }
        }
    }

    recombPos.resize(dimn1, std::vector<std::vector<int>>(dimn2, std::vector<int>(1, 0)));
}

MI_matrix::~MI_matrix()
{
    for (int i = 0; i < dimn1; i++)
    {
        for (int j = 0; j < dimn2; j++)
        {
            delete[] matrixMI[i][j];
        }
        delete[] matrixMI[i];
    }
    delete[] matrixMI;
}

std::vector<int> MI_matrix::getMI(int cur_iter_n, int cur_meiosis_n)
{

    std::vector<int> result;

    result.push_back(matrixMI[cur_iter_n - 1][cur_meiosis_n - 1][0]);
    result.push_back(matrixMI[cur_iter_n - 1][cur_meiosis_n - 1][1]);
    result.push_back(matrixMI[cur_iter_n - 1][cur_meiosis_n - 1][2]);

    return result;
}

void MI_matrix::setRecombPos(int cur_iter_n, int cur_meiosis_n, int num_switches, std::vector<int> recomb_pos)
{

    if (num_switches > 0)
    {
        recombPos[cur_iter_n - 1][cur_meiosis_n - 1].erase(recombPos[cur_iter_n - 1][cur_meiosis_n - 1].begin());

        recombPos.push_back(vector<vector<int>>());
        recombPos[cur_iter_n - 1].push_back(vector<int>());

        for (int i : recomb_pos)
            recombPos[cur_iter_n - 1][cur_meiosis_n - 1].push_back(i);
    }
}

void MI_matrix::printRecombPos(int cur_iter_n, int cur_meiosis_n, int elem_num)
{
    std::cout << recombPos[cur_iter_n - 1][cur_meiosis_n - 1][elem_num] << std::endl;
}

std::vector<int> MI_matrix::getRecombPos(int cur_iter_n, int cur_meiosis_n)
{

    std::vector<int> recomb_positions;

    for (std::vector<int>::iterator iter1 = recombPos[cur_iter_n - 1][cur_meiosis_n - 1].begin();
         iter1 != recombPos[cur_iter_n - 1][cur_meiosis_n - 1].end(); iter1++)
    {
        recomb_positions.push_back(*iter1);
    }

    if (!recomb_positions.empty())
    {
        return recomb_positions;
    }

    else

    {
        recomb_positions.push_back(0);
        return recomb_positions;
    }
}

void MI_matrix::setVal(int cur_iter_n, int cur_meiosis_n, int chr_origin, int start_chr, int num_switches)
{
    matrixMI[cur_iter_n - 1][cur_meiosis_n - 1][0] = chr_origin;
    matrixMI[cur_iter_n - 1][cur_meiosis_n - 1][1] = start_chr;
    matrixMI[cur_iter_n - 1][cur_meiosis_n - 1][2] = num_switches;
}

void MI_matrix::printMeiosisLine(int cur_iter_n, int cur_meiosis_n)
{
    std::cout << cur_iter_n << "  " << cur_meiosis_n << "  " << matrixMI[cur_iter_n - 1][cur_meiosis_n - 1][0]
         << "  " << matrixMI[cur_iter_n - 1][cur_meiosis_n - 1][1] << "  " << matrixMI[cur_iter_n - 1][cur_meiosis_n - 1][2] << std::endl;
}


class FrameworkMarkerPos
{

public:
    double getPos_cM(int);

    void setVal(double);

    std::tuple<char, double, double, int, int> getFlanking_Markers_Pos_cM(double);

    int getMarker_Num(double);

    std::vector<double> getFrameMarkers4ROI(std::vector<double> &);

    FrameworkMarkerPos() : framework_marker_pos() {}

private:
    std::vector<double> framework_marker_pos;
    std::vector<double>::iterator iter1;
};

double FrameworkMarkerPos::getPos_cM(int marker_num)
{

    return framework_marker_pos[marker_num - 1];
}

void FrameworkMarkerPos::setVal(double marker_cM_pos)
{

    framework_marker_pos.push_back(marker_cM_pos);
}

std::tuple<char, double, double, int, int> FrameworkMarkerPos::getFlanking_Markers_Pos_cM(double marker_cM_pos)
{

    iter1 = std::find(framework_marker_pos.begin(), framework_marker_pos.end(), marker_cM_pos);
    if (iter1 != framework_marker_pos.end())
    {
        return std::make_tuple('E', *iter1, *iter1, getMarker_Num(*iter1), getMarker_Num(*iter1));
    }

    else if (framework_marker_pos[0] > marker_cM_pos)
    {
        return std::make_tuple('N', -10, framework_marker_pos[0], 0,
        getMarker_Num(framework_marker_pos[0]));
    }

    else if (framework_marker_pos[framework_marker_pos.size() - 1] < marker_cM_pos)
    {
        return std::make_tuple('N', framework_marker_pos[framework_marker_pos.size() - 1], -10,
        getMarker_Num(framework_marker_pos[framework_marker_pos.size() - 1]), 0);
    }

    else
    {
        for (std::vector<double>::size_type i = 0; i < framework_marker_pos.size() - 1; i++)
        {
            if (framework_marker_pos[i] < marker_cM_pos &&
                framework_marker_pos[i + 1] > marker_cM_pos)
            {
                return std::make_tuple('N', framework_marker_pos[i], framework_marker_pos[i + 1],
                getMarker_Num(framework_marker_pos[i]), getMarker_Num(framework_marker_pos[i + 1]));
            }
        }
    }

    return std::make_tuple ('A', 0, 0, 0, 0);
}

int FrameworkMarkerPos::getMarker_Num(double marker_cM_pos)
{

    int marker_num = -1;
    iter1 = std::find(framework_marker_pos.begin(), framework_marker_pos.end(), marker_cM_pos);
    if (iter1 != framework_marker_pos.end())
    {
        marker_num = std::distance(framework_marker_pos.begin(), iter1) + 1;
        return marker_num;
    }
    return -1;
}

std::vector<double> FrameworkMarkerPos::getFrameMarkers4ROI(std::vector<double> &linkage_rg)
{

    std::vector<double> result;

    for (std::vector<double>::iterator iter1 = framework_marker_pos.begin();
         iter1 != framework_marker_pos.end(); iter1++)
    {

        if (*iter1 >= linkage_rg[0] && *iter1 <= linkage_rg[1])
        {
            result.push_back(*iter1);
        }
    }

    return result;
}


class DenseMarkerGenos
{

public:
    DenseMarkerGenos();

    void setGeno(int, int, std::vector<std::string>);

    std::vector<std::string> getGeno(int, int);

private:
    DenseGenoMap dense_genotypes;

};

DenseMarkerGenos::DenseMarkerGenos() : dense_genotypes{}
{
}

void DenseMarkerGenos::setGeno(int var_bp_pos, int subject_index, std::vector<std::string> al_1_2)
{

    dense_genotypes[var_bp_pos][subject_index].push_back(al_1_2[0]);
    dense_genotypes[var_bp_pos][subject_index].push_back(al_1_2[1]);

}

std::vector<std::string> DenseMarkerGenos::getGeno(int var_bp_pos, int subject_index)
{
    std::vector<std::string> genotype;

    std::map<int, std::map<int, std::vector<std::string>>>::iterator iter1;
    std::map<int, std::vector<std::string>>::iterator iter2;

    iter1 = dense_genotypes.find(var_bp_pos);
    iter2 = dense_genotypes.at((*iter1).first).find(subject_index);

    std::vector<std::string> temp1 = (*iter2).second;
   
    genotype.push_back(temp1[0]);
    genotype.push_back(temp1[1]);

    return genotype;
}

double calcWeight(int m1, int m2, double d1, double d2, double pos)
{
    double p = 0.0;
    double x = (pos - d1) / 100;
    double y = (d2 - pos) / 100;
    double r_x = (1 - exp(-2 * x)) / 2;
    double r_y = (1 - exp(-2 * y)) / 2;
    double r_xy = (1 - exp(-2 * (x + y))) / 2;

    if ((d1 <= pos) && (pos < d2))
    {
        if ((m1 == 0) && (m2 == 0))
        {
            p = r_x * r_y / (1 - r_xy);
        }
        else if ((m1 == 0) && (m2 == 1))
        {
            p = r_x * (1 - r_y) / (r_xy);
        }
        else if ((m1 == 1) && (m2 == 0))
        {
            p = (1 - r_x) * r_y / (r_xy);
        }
        else if ((m1 == 1) && (m2 == 1))
        {
            p = (1 - r_x) * (1 - r_y) / (1 - r_xy);
        }
    }

    else if (d2 == -10)
    {

        if (m1 == 0)
        {
            p = r_x;
        }
        else if (m1 == 1)
        {
            p = (1 - r_x);
        }
    }

    else if (d1 == -10)
    {
        if (m2 == 0)
        {
            p = r_y;
        }
        else if (m2 == 1)
        {
            p = (1 - r_y);
        }
    }
    return p;
}


class FGL_matrix
{
public:
    FGL_matrix(int, int, Pedigree_info *, MI_matrix *, FrameworkMarkerPos *,
    std::mt19937&, DenseMarkerPos *);

    ~FGL_matrix();

    void setFGL(int);

    void setFGLfm(double);

    std::vector<int> getFGL(int, int);

    int estimateMI(double, int, int);

    int getMI_FrameMarker(int, int, int);

private:  
    int dimn1, dimn2;  
    Pedigree_info *ped1 = nullptr;
    MI_matrix *mi_matrix = nullptr;
    FrameworkMarkerPos *frame_m_pos = nullptr;
    std::mt19937& rng;
    DenseMarkerPos *dense_m_pos = nullptr;
    int ***matrixFGL = nullptr;
};

FGL_matrix::FGL_matrix(int num_iter, int num_subjects,
Pedigree_info *cur_ped, MI_matrix *mi_matrix_ptr,
FrameworkMarkerPos *framework_marker_pos, std::mt19937& rngIn,
DenseMarkerPos *dense_marker_pos) : dimn1{num_iter},
dimn2{num_subjects}, ped1{cur_ped}, mi_matrix{mi_matrix_ptr},
frame_m_pos{framework_marker_pos}, rng{rngIn}, dense_m_pos{dense_marker_pos}
{

    matrixFGL = new int **[dimn1];
    for (int i = 0; i < dimn1; i++)
    {
        matrixFGL[i] = new int *[dimn2];
        for (int j = 0; j < dimn2; j++)
        {
            matrixFGL[i][j] = new int[2];
            for (int k = 0; k < 2; k++)
            {
                matrixFGL[i][j][k] = -1;
            }
        }
    }
}

FGL_matrix::~FGL_matrix()
{
    for (int i = 0; i < dimn1; i++)
    {
        for (int j = 0; j < dimn2; j++)
        {
            delete[] matrixFGL[i][j];
        }
        delete[] matrixFGL[i];
    }
    delete[] matrixFGL;
}

void FGL_matrix::setFGL(int variant_bp)
{

    double dense_m_cM = dense_m_pos->getPos_cM(variant_bp);

    for (int k = 0; k < dimn1; k++)
    {

        for (size_t i = 0; i < ped1->cur_ped.size(); i++)
        {

            if (ped1->isFounder(i) == 1)
            {

                matrixFGL[k][i][0] = ped1->cur_ped[i].mat_fgl;
                matrixFGL[k][i][1] = ped1->cur_ped[i].pat_fgl;
            }
            else
            {
                if (ped1->cur_ped[i].mat_mi != 0 && ped1->cur_ped[i].pat_mi != 0)
                {
                    int m_ind = ped1->getSubjectIndex(ped1->cur_ped[i].mother_id);
                    int p_ind = ped1->getSubjectIndex(ped1->cur_ped[i].father_id);

                    int m_mi = estimateMI(dense_m_cM, k + 1, ped1->cur_ped[i].mat_mi);
                    int p_mi = estimateMI(dense_m_cM, k + 1, ped1->cur_ped[i].pat_mi);

                    
                    int m_fgl = matrixFGL[k][m_ind][m_mi];
                    int p_fgl = matrixFGL[k][p_ind][p_mi];

                    matrixFGL[k][i][0] = m_fgl;
                    matrixFGL[k][i][1] = p_fgl;
                }
                else
                {
                    if (ped1->cur_ped[i].mat_mi == 0)
                    {

                        int m_ind = ped1->getSubjectIndex(ped1->cur_ped[i].mother_id);
                        int p_ind = ped1->getSubjectIndex(ped1->cur_ped[i].father_id);

                        int p_mi = estimateMI(dense_m_cM, k + 1, ped1->cur_ped[i].pat_mi);

                        int m_fgl = matrixFGL[k][m_ind][0];
                        int p_fgl = matrixFGL[k][p_ind][p_mi];

                        matrixFGL[k][i][0] = m_fgl;
                        matrixFGL[k][i][1] = p_fgl;
                    }
                    else
                    {

                        int m_ind = ped1->getSubjectIndex(ped1->cur_ped[i].mother_id);
                        int p_ind = ped1->getSubjectIndex(ped1->cur_ped[i].father_id);

                        int m_mi = estimateMI(dense_m_cM, k + 1, ped1->cur_ped[i].mat_mi);


                        int m_fgl = matrixFGL[k][m_ind][m_mi];
                        int p_fgl = matrixFGL[k][p_ind][0];

                        matrixFGL[k][i][0] = m_fgl;
                        matrixFGL[k][i][1] = p_fgl;
                    }
                }
            }

        }
    }
}

void FGL_matrix::setFGLfm(double dense_m_cM)
{

    for (int k = 0; k < dimn1; k++)
    {

        for (size_t i = 0; i < ped1->cur_ped.size(); i++)
        {

            if (ped1->isFounder(i) == 1)
            {

                matrixFGL[k][i][0] = ped1->cur_ped[i].mat_fgl;
                matrixFGL[k][i][1] = ped1->cur_ped[i].pat_fgl;
            }
            else
            {
                if (ped1->cur_ped[i].mat_mi != 0 && ped1->cur_ped[i].pat_mi != 0)
                {
                    int m_ind = ped1->getSubjectIndex(ped1->cur_ped[i].mother_id);
                    int p_ind = ped1->getSubjectIndex(ped1->cur_ped[i].father_id);

                    int m_mi = estimateMI(dense_m_cM, k + 1, ped1->cur_ped[i].mat_mi);
                    int p_mi = estimateMI(dense_m_cM, k + 1, ped1->cur_ped[i].pat_mi);

                    int m_fgl = matrixFGL[k][m_ind][m_mi];
                    int p_fgl = matrixFGL[k][p_ind][p_mi];

                    matrixFGL[k][i][0] = m_fgl;
                    matrixFGL[k][i][1] = p_fgl;
                }
                else
                {
                    if (ped1->cur_ped[i].mat_mi == 0)
                    {

                        int m_ind = ped1->getSubjectIndex(ped1->cur_ped[i].mother_id);
                        int p_ind = ped1->getSubjectIndex(ped1->cur_ped[i].father_id);

                        int p_mi = estimateMI(dense_m_cM, k + 1, ped1->cur_ped[i].pat_mi);

                        int m_fgl = matrixFGL[k][m_ind][0];
                        int p_fgl = matrixFGL[k][p_ind][p_mi];

                        matrixFGL[k][i][0] = m_fgl;
                        matrixFGL[k][i][1] = p_fgl;
                    }
                    else
                    {

                        int m_ind = ped1->getSubjectIndex(ped1->cur_ped[i].mother_id);
                        int p_ind = ped1->getSubjectIndex(ped1->cur_ped[i].father_id);

                        int m_mi = estimateMI(dense_m_cM, k + 1, ped1->cur_ped[i].mat_mi);

                        int m_fgl = matrixFGL[k][m_ind][m_mi];
                        int p_fgl = matrixFGL[k][p_ind][0];

                        matrixFGL[k][i][0] = m_fgl;
                        matrixFGL[k][i][1] = p_fgl;
                    }
                }
            }
        }
    }
}

std::vector<int> FGL_matrix::getFGL(int iter_num, int subject_index)
{

    std::vector<int> result;

    result.push_back(matrixFGL[iter_num - 1][subject_index][0]);
    result.push_back(matrixFGL[iter_num - 1][subject_index][1]);

    return result;
}

int FGL_matrix::estimateMI(double dense_m_cM, int cur_iter_num,
int MI_line_num)
{

    char frame_m_config;
    double left_fl_m_cM;
    double right_fl_m_cM;
    int left_fl_m_num;
    int right_fl_m_num;

    tie(frame_m_config, left_fl_m_cM, right_fl_m_cM, left_fl_m_num, right_fl_m_num) =
    frame_m_pos->getFlanking_Markers_Pos_cM(dense_m_cM);

    if (frame_m_config == 'N')
    {

        int frame_m_mi_left = getMI_FrameMarker(cur_iter_num, MI_line_num, left_fl_m_num);
        int frame_m_mi_right = getMI_FrameMarker(cur_iter_num, MI_line_num, right_fl_m_num);

        double prob_inh_father = calcWeight(frame_m_mi_left, frame_m_mi_right, left_fl_m_cM,
        right_fl_m_cM, dense_m_cM);

        double randNum = rand(rng);
        int dense_m_mi;

        if (randNum <= prob_inh_father)
        {
            dense_m_mi = 1;
        }
        else
        {
            dense_m_mi = 0;
        }

        return dense_m_mi;
    }

    else if (frame_m_config == 'E')
    {
        int dense_m_mi = getMI_FrameMarker(cur_iter_num, MI_line_num, left_fl_m_num);
        return dense_m_mi;
    }

    return -1;
}

int FGL_matrix::getMI_FrameMarker(int cur_iter_num, int MI_line_num, int frame_m_num)
{

    int result = -1;
    std::vector<int> recomb_pos;

    std::vector<int> mi_values;
    mi_values = mi_matrix->getMI(cur_iter_num, MI_line_num);

    if (mi_values[2] == 0)
    {
        result = mi_values[1];
    }

    else

    {
        recomb_pos = mi_matrix->getRecombPos(cur_iter_num, MI_line_num);

        if (frame_m_num < recomb_pos[recomb_pos.size() - 1] && frame_m_num >= recomb_pos[0])
        {

            for (std::vector<double>::size_type i = 0; i < recomb_pos.size() - 1; i++)
            {

                if (frame_m_num >= recomb_pos[i] && frame_m_num < recomb_pos[i + 1])
                {

                    int cur_switch_num = i + 1;

                    if (cur_switch_num % 2 == 0)
                    {
                        result = mi_values[1];
                    }
                    else
                    {
                        if (mi_values[1] == 0)
                        {
                            result = 1;
                        }

                        else
                        {
                            result = 0;
                        }
                    }
                    goto LB1;
                }
            }
        }

        else if (frame_m_num >= recomb_pos[recomb_pos.size() - 1])
        {

            int cur_switch_num = recomb_pos.size();

            if (cur_switch_num % 2 == 0)
            {
                result = mi_values[1];
            }
            else
            {
                if (mi_values[1] == 0)
                {
                    result = 1;
                }

                else
                {
                    result = 0;
                }
            }
            goto LB1;
        }

        else if (frame_m_num < recomb_pos[0])
        {
            result = mi_values[1];
        }
    }

LB1:;

    return result;
}


std::vector<string> split2(const std::string &s, char delim)
{
    std::vector<string> elems;
    std::stringstream ss(s);
    std::string number;
    while (std::getline(ss, number, delim))
    {
        elems.push_back(number);
    }
    return elems;
}


class FGLshareGroup2
{

public:
    
    FGLshareGroup2(Pedigree_info *, Parameters *, MI_matrix *,
    FrameworkMarkerPos *, DenseMarkerPos *, std::mt19937&);
    
    std::vector<int> evalFGLshareMain();

    void getROIFrameMarkers();

    void compFGLshare(int, int, std::vector<std::vector<std::vector<int>>> &);

    std::vector<std::string> keepCases(std::vector<std::string> &);
    
    void recordFGLshareIter(int, int, std::vector<string> &, std::vector<std::vector<std::vector<int>>> &);

    void selectShareGrp4FGL(std::vector<std::vector<std::vector<int>>> &, std::map<int, std::vector<int>> &, std::map<int, std::vector<int>> &);

    std::vector<int> detectContShareRegion(std::vector<int> &);

    int detIdShareEl (std::vector<int>, std::vector<int>);

    void evalFGLshareLvF(int, std::map<int, std::vector<int>>&, std::map<int, std::vector<int>>&, std::map<std::string, std::vector<int>>&);
    
    std::vector<int> evalFGLshareLv3(std::map<std::string, std::vector<int>> &);

private:
    std::vector<double> roiMarkers; // cM positions for linkage markers in ROI
    std::vector<double> linkage_boundaries;
    Pedigree_info *ped1 = nullptr;
    Parameters *pars = nullptr;
    MI_matrix *mi_matrix = nullptr;
    std::mt19937& rng;
    FrameworkMarkerPos *frame_m_pos = nullptr;
    DenseMarkerPos *dense_m_pos = nullptr;

    int mapKeyNum;
    int mapKeyNum2;
    int mapKeyNum3;
    int num_chrs, num_iters;
    int minRg;// min sharing region equals 25% of ROI
    double maxLODmPoscM;// cM position of linkage marker with maxLOD score
    int maxLODmPos;// sequential position of linkage marker with LODmax score in ROI

    std::vector<std::vector<int>> equilClasses; // list of iteration numbers where sharing for same group of cases was detected

    std::vector<int> cases;
};

FGLshareGroup2::FGLshareGroup2(Pedigree_info *current_ped_ptr,
Parameters *current_pars_ptr, MI_matrix *mi_matrix_ptr,
FrameworkMarkerPos *framework_marker_pos, DenseMarkerPos *dense_marker_pos, std::mt19937& rngIn) : 
ped1{current_ped_ptr}, pars{current_pars_ptr}, mi_matrix{mi_matrix_ptr}, rng{rngIn},
frame_m_pos{framework_marker_pos}, dense_m_pos{dense_marker_pos}
{

    linkage_boundaries = pars->linkage_region;
    num_iters = pars->num_iter;
    num_chrs = ped1->getNumSubjects() * 2;
    mapKeyNum = -1;
    mapKeyNum2 = -1; // key for equilClasses structure
    mapKeyNum3 = -1;
    maxLODmPoscM = pars->maxLODmarker;
    
}

/// main function for this class
std::vector<int> FGLshareGroup2::evalFGLshareMain()
{

    this->getROIFrameMarkers();

    minRg = round ((roiMarkers.size()) * pars->minShareLength); // minimum length of shared region in number of linkage markers, default is 25% of the number of linkage markers in ROI
  

    // identifies vector element # containing value equal to cM position of maxLOD linkage marker
    for (size_t i =1; i < roiMarkers.size()-1; i++)
    {
        if (maxLODmPoscM > roiMarkers[i-1] && maxLODmPoscM < roiMarkers[i+1])
        {
            maxLODmPos = i;
        }
    }


    std::map<std::string, std::vector<int>> sharingGrpsTotal; // key: combination of ordered case IDs sharing same FGL; values: [0] = mapKeyNum2, [1]+ = IDs of cases sharing same FGL

    for ( int j=1; j<=num_iters; j++) // real iteration numbers starting with 1
    {

        std::vector<std::vector<std::vector<int>>> sharing4FGL; // v1 - fgl # (real fgl-1); v2 - linkage maker # (marker #-1) in ROI; v3 - some info (0= # of sharing cases with WGS, 1= # of all sharing cases, 2+ = all sharing cases' ids)
        int num_fgls = ped1->getNumFounders() * 2;
        sharing4FGL.resize(num_fgls);
        for (std::vector<std::vector<std::vector<int>>>::iterator it = sharing4FGL.begin(); it != sharing4FGL.end(); it++)
        {
            it->resize(roiMarkers.size());
        }

        std::map<int, std::vector<int>> sharingGrpsV2; // keys: arbitrary sequential numbers; values: el0 = fgl; el1 = # of sharing cases with WGS; el2 = # of total sharing cases; el3 + = sharing cases' ids

        std::map<int, std::vector<int>> sharingGrpsRgV2; // keys: arbitrary sequential numbers; values: el0 = # of markers within continues sharing region; el1 += linkage marker positions in ROI where sharing has occured

        for (size_t i = 0; i < roiMarkers.size(); i++)
        {
            this->compFGLshare(i, j, sharing4FGL);
        } // end of ROI linkage markers loop

        this->selectShareGrp4FGL(sharing4FGL, sharingGrpsV2, sharingGrpsRgV2);

        this->evalFGLshareLvF(j, sharingGrpsV2, sharingGrpsRgV2, sharingGrpsTotal);

        mapKeyNum = -1;
        mapKeyNum3 = -1;
       
    } // end of iteration loop


    std::vector<int> result1 = this->evalFGLshareLv3(sharingGrpsTotal); // el0 = iter# for phasing; el1 += list of cases from the core set

    return result1;
}

/// gets cM positions for frame markers in ROI
void FGLshareGroup2::getROIFrameMarkers()
{
    roiMarkers = frame_m_pos->getFrameMarkers4ROI(linkage_boundaries);
}

/// identifies groups of cases sharing the same FGL at a linkage marker position in ROI
void FGLshareGroup2::compFGLshare(int m_el, int iter_num,
std::vector<std::vector<std::vector<int>>> &sharing4FGL)
{

    int num_fgls = ped1->getNumFounders() * 2;

    FGL_matrix fgl_matrix(num_iters, ped1->getNumSubjects(),
    ped1, mi_matrix, frame_m_pos, rng, dense_m_pos);// initiates FGL matrix for a linkage marker acroos all iteration

    fgl_matrix.setFGLfm(roiMarkers[m_el]); // generates FGL matrix spanning all iterations for a linkage marker based on its cM_pos

    FGL_matrix *fgl_matrix_ptr = &fgl_matrix;

    for (int fgl_val = 0; fgl_val < num_fgls; fgl_val++) // loop over all FGL values
    {

        std::vector<std::string> temp_same_fgl_subjects;

        for (int sub_ind = 0; sub_ind < ped1->getNumSubjects(); sub_ind++)
        {

            std::vector<int> fgls_cur = fgl_matrix_ptr->getFGL(iter_num, sub_ind);
            
            if (fgl_val == fgls_cur[0] || fgl_val == fgls_cur[1])
            {

                std::string sub_id = ped1->cur_ped[sub_ind].subject_id;
                temp_same_fgl_subjects.push_back(sub_id);
            }

        } // sub_ind loop

        if (temp_same_fgl_subjects.size() > 1)
        {

            std::vector<std::string> same_fgl_cases = this->keepCases(temp_same_fgl_subjects);// removes ids who are not cases

            if (same_fgl_cases.size() > 1)
            {
                this->recordFGLshareIter(m_el, fgl_val, same_fgl_cases, sharing4FGL);
            }
        }

    } // fgl_val loop
}

/// selects only ids which are cases
std::vector<std::string> FGLshareGroup2::keepCases(std::vector<std::string> &subjects_all)
{

    std::vector<std::string> results;

    for (size_t i = 0; i < subjects_all.size(); i++)
    {

        int SubInd = ped1->getSubjectIndex(subjects_all[i]);

        if (ped1->cur_ped[SubInd].pheno == 2)
        {
            results.push_back(subjects_all[i]);
        }
    }

    std::sort(results.begin(), results.end());

    return results;// returns only cases' ids
}

/// stores FGL sharing information at a linkage marker in ROI
void FGLshareGroup2::recordFGLshareIter(int m_el, int fgl, std::vector<string> &same_fgl_cases, std::vector<std::vector<std::vector<int>>> &sharing4FGL)
{

    std::vector<int> temp3;
    temp3.push_back(ped1->getNumWGSlist(same_fgl_cases));// # of sharing cases with WGS data
    temp3.push_back(same_fgl_cases.size());// total # of sharing cases

    for (size_t i = 0; i < same_fgl_cases.size(); i++)
    {
        temp3.push_back(stoi(same_fgl_cases[i]));// ids of all sharing cases
    }

    int fgl_el = fgl - 1;

    for (size_t i = 0; i < temp3.size(); i++)
    {
        sharing4FGL[fgl_el][m_el].push_back(temp3[i]);
    }
}

/// identifies largest group of subjects sharing same FGL along longest contiguous region in ROI
void FGLshareGroup2::selectShareGrp4FGL(std::vector<std::vector<std::vector<int>>> &sharing4FGL, std::map<int, std::vector<int>> &sharingGrpsV2, std::map<int, std::vector<int>> &sharingGrpsRgV2)
{

    for (size_t i = 0; i < sharing4FGL.size(); i++) // fgl
    {

        std::map<int, std::vector<int>> fglShareIds; // key = sequential number; vector - sharing ids
        std::map<int, std::vector<int>> sharePos;    // key = sequential number; vector - positions for multiple sharing regions in ROI starting from 0
        std::map<int, std::vector<int>> sharePosFinal; // key = sequential number; vector - positions for largest contiguous sharing region
        
        int map_key = 0;// key for fglShareIds, sharePos, and sharePosFinal maps

        std::vector <int> ids2select;// contains # of cases for each unique grp sharing FGL
        std::vector <int> shareRg2select;// contains # of markers within largest contiguous sharing region shared by subjects whose number is in ids2select vector with same element #
        std::vector <int> keyHolder;// contains map keys associated with each el in dis2select and shareRg2select vecotrs

        for (size_t j = 0; j < sharing4FGL[i].size(); j++) // sequential # of a linkage marker in ROI
        {
            
            std::vector<int> shareGrpInit; // subjects sharing ith FGL at j linkage marker position

            if ((sharing4FGL[i][j].size() - 2)> 1) // -2 because first 2 elements of vector contain different type of information; this insures that we have at least two people sharing same FGL
            {
                for (size_t k = 2; k < sharing4FGL[i][j].size(); k++) // starting from 2 because 0 - # of sharing cases with WGS; 1 - # of total sharing cases (with and without WGS)
                {
                    shareGrpInit.push_back(sharing4FGL[i][j][k]);
                }
            }

                for (size_t n = 0; n < sharing4FGL[i].size(); n++) // sequential # of a linkage marker in ROI
                {

                    std::vector<int> sharePosCur;// subjects sharing ith FGL at nth linkage marker position

                    if ((sharing4FGL[i][n].size() -2) > 1)
                    {
                        for (size_t k = 2; k < sharing4FGL[i][n].size(); k++) // starting from 2 because 0 - cased with WGS; 1 - cases total (with and without WGS)
                        {
                            sharePosCur.push_back(sharing4FGL[i][n][k]);
                        }
                    }

                        std::vector<int> sameFGLshare;
                        if (!shareGrpInit.empty() && !sharePosCur.empty())
                        {
                            Comp_arrays_venn<int> result1(&shareGrpInit, &sharePosCur);
                            sameFGLshare = result1.getIntersection();
                        }

                        
                        ///////////// record result
                        if (!sameFGLshare.empty())
                        {

                            for (std::map<int, std::vector<int>>::iterator iter1 = fglShareIds.begin(); iter1 != fglShareIds.end(); iter1++)
                            {

                                std::vector<int> temp1 = iter1->second;

                                if (temp1 == sameFGLshare)
                                {
                                    int cur_key = iter1->first;
                                    sharePos[cur_key].push_back(n);
                                    goto LB10;
                                }
                            }

                            map_key++;// starts over for each fgl from 1
                            fglShareIds.insert(make_pair(map_key, sameFGLshare));
                            sharePos[map_key].push_back(n);
                        }
                    LB10:;
                }
            
        } // end of ROI linkage markers loop


        // remove duplicate linkage marker positions where sharing for a group of subjects has occured
        for (std::map<int, std::vector<int>>::iterator iter1 = sharePos.begin(); iter1 != sharePos.end(); iter1++)
        {

            std::vector<int> tmp;// contains linkage marker positions where sharing for a group of subjects has occured
            int cur_key = iter1->first;
            tmp = iter1->second;

            sort(tmp.begin(), tmp.end());

            vector<int>::iterator tmp_iter;
            tmp_iter = std::unique(tmp.begin(), tmp.end());
            tmp.resize(std::distance(tmp.begin(), tmp_iter));

            // check if sharing contiguous

            std::vector <int> res1 =  detectContShareRegion (tmp);// returns linkage marker positions for the longest contiguous sharing or -1 if no sharing detected


                if (fglShareIds[cur_key].size() > 1 && res1[0] != -1)
                {

                    std::vector<int> check = fglShareIds[cur_key];

                    keyHolder.push_back (cur_key);
                    ids2select.push_back(check.size());

                    std::vector <int> temp3;
                    temp3.push_back(res1.size());

                    for (size_t i = 0; i < res1.size(); i++)
                    {
                        temp3.push_back(res1[i]);         
                    }
                    sharePosFinal.insert(make_pair(cur_key, temp3));

                    shareRg2select.push_back(res1.size());
                        
                      
                }

        }



        if (ids2select.size() >= 2)
        {
            
            int el_num = detIdShareEl (ids2select, shareRg2select); // 1st vector - size of FGL sharing group; 2nd vector - length of FGL sharing region in ROI

            if (el_num == -1)
            {
                goto LB2;
            }

        /////// record selected sharing group


            mapKeyNum3++;
            int map_key2 = mapKeyNum3 + 1;
            std::vector<int> temp1;

            std::vector<std::string> tp1;
            std::vector<int> shareGrp = fglShareIds[keyHolder[el_num]];
            for (size_t i = 0; i < shareGrp.size(); i++)
            {
                tp1.push_back(std::to_string(shareGrp[i]));
            }

            temp1.push_back(i + 1); // push FGL #
            temp1.push_back(ped1->getNumWGSlist(tp1)); // # of subjects with WGS data sharing same FGL
            temp1.push_back(shareGrp.size()); // total # of subjects sharing same FGL

            for (size_t i = 0; i < shareGrp.size(); i++)
            {
                temp1.push_back(shareGrp[i]); // subjects sharing same FGL
            }
            sharingGrpsV2.insert(make_pair(map_key2, temp1));
                
        //// record unique contiguous sharing
            std::vector<int> sharePs = sharePosFinal[keyHolder[el_num]];
            sharingGrpsRgV2.insert(make_pair(map_key2, sharePs));

        }
        LB2:;

        // if only one record of FGL sharing in iteration
        if (ids2select.size() == 1 && shareRg2select[0] >= minRg)
        {

                    int el_num = 0;
                    mapKeyNum3++;
                    int map_key2 = mapKeyNum3 + 1;
                    std::vector<int> temp1;

                    std::vector<std::string> tp1; // subjects sharing same FGL
                    std::vector<int> shareGrp = fglShareIds[keyHolder[el_num]];
                    for (size_t i = 0; i < shareGrp.size(); i++)
                    {
                        tp1.push_back(std::to_string(shareGrp[i]));
                    }

                    temp1.push_back(i + 1); // push FGL #
                    temp1.push_back(ped1->getNumWGSlist(tp1)); // # of subjects in shared group with WGS data
                    temp1.push_back(shareGrp.size()); // total # of subjects in shared group

                    for (size_t i = 0; i < shareGrp.size(); i++)
                    {
                        temp1.push_back(shareGrp[i]); // sujbects sharing same FGL
                    }
                    sharingGrpsV2.insert(make_pair(map_key2, temp1)); // map key is +1 from el of vector
                

            // record unique contiguous sharing
            std::vector<int> sharePs = sharePosFinal[keyHolder[el_num]];
            sharingGrpsRgV2.insert(make_pair(map_key2, sharePs));

        }

    } // end of fgl loop
      
}

/// determines a region with contiguous sharing in ROIss
std::vector<int> FGLshareGroup2::detectContShareRegion(std::vector<int> &sharePos)
{
// determines longest contiguous sharing region
// returns vector with marker positions at which sharing has occured or -1 if no sharing determined
    std::vector<int> result1; // if el 0 = -1 then negative result
    std::sort(sharePos.begin(), sharePos.end()); // linkage marker numbers start from 0
    std::vector<int> share_bounds; // contains element numbers for vector sharePos to define boundaries of contiguous shared regions in ROI
    std::vector<int> share_size; // contains length of contiguous shared regions in numbers of linkage markers
    std::vector<int> selectedShare; // contains elements (two for each region) for sharePos vector to define boundaries (start and end) of shared region
    int maxLength;
    int status1 = 0;

    if (sharePos.size() < 2)
    {
        result1.push_back(-1);
        goto LB1;
    }

    // determines start and end of contiguous shared regions in ROI; should contain two entries per shared region
    for (size_t i = 0; i < sharePos.size()-1; i++)
    {

        if ((sharePos[i] == (sharePos[i + 1] - 1)) && status1 == 0)
        {
            status1 = 1; // start of sharing region
            share_bounds.push_back(i);
        }

        else if ((sharePos[i] == (sharePos[i + 1] - 1)) && status1 == 1)
        {
           // moving through sharing region
        }

        else if ((sharePos[i] != (sharePos[i + 1] - 1)) && status1 == 1)
        {
            // end of sharing region
            share_bounds.push_back(i);
            status1 = 0;

        }

        else if ((sharePos[i] != (sharePos[i + 1] - 1)) && status1 == 0)
        {
            // nothing, keep going
        }

    }

    if (status1 == 1)
    {
        share_bounds.push_back(sharePos.size()-1);
        status1 = 0;
    }

    if (share_bounds.empty())
    {
        result1.push_back(-1);
        goto LB1;
    }


    // determines length of contiguous shared region
    for (size_t i = 1; i < share_bounds.size(); i++)
    {
        share_size.push_back(share_bounds[i] - share_bounds[i - 1]);
    }

    
    if (share_size.size() > 1)
    {
        maxLength = *max_element(share_size.begin(), share_size.end());
    }
    else
    {
        maxLength = share_size[0];
    }


    for (size_t i = 0; i < share_size.size(); i++)
    {

        if (share_size[i] == maxLength)
        {
            selectedShare.push_back(share_bounds[2*i]);
            selectedShare.push_back(share_bounds[(2*i) + 1]);
        }
    }

    // outputs linkage marker numbers from largest shared region in ROI into results vector
    for (size_t i = 0; i < sharePos.size(); i++)
    {
        if (selectedShare[0] >= 0 && selectedShare[1] >= 0 &&
            i >= static_cast<size_t>(selectedShare[0]) &&
            i <= static_cast<size_t>(selectedShare[1]))
        {
            result1.push_back(sharePos[i]);
        }
    }

    // checks whether selected region contains linkage marker with maxLOD score
    for (size_t i =0; i< result1.size(); i++)
    {
      if (result1[i] == maxLODmPos)
      {
        goto LB1;
      }
    }

    // if selected sharing region does not contain linkage marker with maxLOD score then negative result
    result1.clear();
    result1.push_back(-1);


    LB1:;
    return result1; // if el0 = -1 then negative result
}

/// determines longest sharing region in ROI and largest sharing group among different options for an FGL
int FGLshareGroup2::detIdShareEl (std::vector<int> ids2select, std::vector<int> shareReg2select) 
{
// this function selects elements in arrays based on min sharing length (>=) computed as 25% of ROI region and then element
// with max number of sharing cases; if more than 1 with max # of cases sharing same fgl
// then pick the one with longest sharing region
// returns element number for arrays
// !!! check whether vector sizes >= 2 before executing this function

    int result1 = -2;
    int minShare = minRg;
    std::vector<int> maxIdsList;
    std::vector<int> elHolder; // contains element # to connect entries in selectedIds to those in ids2select and shareReg2select

   // selection based on minimum sharing length
    std::vector<int> selectedIds; // nunmber of subjects sharing same FGL selected after passing minimum sharing length threshold
    for (size_t i= 0; i < shareReg2select.size(); i++)
    {
        if (shareReg2select[i] >= minShare)
        {
            selectedIds.push_back(ids2select[i]);
            elHolder.push_back(i);
        }

    }

    // what if none selected based on length of sharing region then negative result and get out
    if (selectedIds.empty())
    {
        result1 = -1;
        goto LB1;
    }

    if (selectedIds.size() == 1)
    {
        /// use one element in selectedIds and get out   
        result1 = elHolder[0];
        goto LB1;
           
    }

    /// further selections based on # of cases in sharing group
    int maxIds;
    if (!selectedIds.empty())
    {
        maxIds = *max_element(selectedIds.begin(), selectedIds.end());
    }


    for (size_t i= 0; i < selectedIds.size(); i++)
    {
        if (selectedIds[i] == maxIds)
        {
           maxIdsList.push_back(elHolder[i]);
        }
    }


    if (maxIdsList.size() == 1)
    {
        result1 = maxIdsList[0];
    }
    else
    { // we get here if there are multiple groups with same maximum # of cases but different sharing length
        std::vector<int> shareLength;

        for (size_t i=0; i< maxIdsList.size(); i++)
        {
            shareLength.push_back(shareReg2select[maxIdsList[i]]);
        }
                
        int maxLength2 = -1;
        if (!shareLength.empty())
        {
            maxLength2 = *max_element(shareLength.begin(), shareLength.end());
        }
            
                
        for (size_t i=0; i< maxIdsList.size(); i++)
        {
            if (shareReg2select[maxIdsList[i]] == maxLength2)
            {
                result1 = maxIdsList[i];
                goto LB1;
            }
        }

    }

LB1:;
  
return result1;
} 


/// selects largest haplotype sharing group of cases with longest haplotype sharing region in ROI for an iteration of IVs
void FGLshareGroup2::evalFGLshareLvF(int iter_num, std::map<int, std::vector<int>> &sharingGrpsV2, std::map<int, std::vector<int>> &sharingGrpsRgV2, std::map<std::string, std::vector<int>> &sharingGrpsTotal)
{

    std::vector<int> keys2keep;
    std::vector<int> keys2keep2;
    std::vector<int> shareLength;
    std::vector<int> casesAll;
    int maxCases = -1;
    int maxLength = -1;

    if (sharingGrpsV2.size() > 1) // selection based on max number of all cases sharing same FGL
    {

        for (std::map<int, std::vector<int>>::iterator iter1 = sharingGrpsV2.begin(); iter1 != sharingGrpsV2.end(); iter1++)
        {
            std::vector<int> temp1 = iter1->second;
            casesAll.push_back(temp1[2]); // push # of all cases
        }

        if (!casesAll.empty())
        {
            maxCases = *max_element(casesAll.begin(), casesAll.end());
        }

        for (std::map<int, std::vector<int>>::iterator iter1 = sharingGrpsV2.begin(); iter1 != sharingGrpsV2.end(); iter1++)
        {
            std::vector<int> temp1 = iter1->second;
            if (temp1[2] == maxCases)
            {
                keys2keep.push_back(iter1->first); // keeps keys for entries with max # of cases sharing same FGL
            }
        }


        if (keys2keep.size() == 1)
        {
        // record result and out

            // this part creates unique map key for group of cases sharing same FGL
            std::string map_key; // key for sharingGrpsTotal structure created by combining sort IDs of cases sharing same FGL

            for (std::map<int, std::vector<int>>::iterator iter1 = sharingGrpsV2.begin(); iter1 != sharingGrpsV2.end(); iter1++)
            {

                int cur_key = iter1->first;

                if (keys2keep[0] == cur_key)
                {

                    std::vector<int> temp1 = iter1->second;

                    for (size_t i = 3; i < temp1.size(); i++)
                    {
                        map_key += std::to_string(temp1[i]);
                    }
                }
            }


            if (sharingGrpsTotal.find(map_key) != sharingGrpsTotal.end())
            {

                for (std::map<std::string, std::vector<int>>::iterator iter1 = sharingGrpsTotal.begin(); iter1 != sharingGrpsTotal.end(); iter1++)
                {
                    if (iter1->first == map_key)
                    {
                        std::vector<int> temp3 = iter1->second;

                        int el2edit = temp3[0]; // this entry equals mapKeyNum2
                        equilClasses[el2edit].push_back(iter_num);

                    }
                }
            }
            else
            {
                mapKeyNum2++;

                std::vector<int> temp4;
                temp4.push_back(mapKeyNum2);

                for (std::map<int, std::vector<int>>::iterator iter1 = sharingGrpsV2.begin(); iter1 != sharingGrpsV2.end(); iter1++)
                {

                    if (iter1->first == keys2keep[0])
                    {

                        std::vector<int> temp1 = iter1->second;
                    
                        for (size_t i = 3; i < temp1.size(); i++)
                        {
                            temp4.push_back(temp1[i]);
                        }
                    
                    }
                }

                sharingGrpsTotal.insert(make_pair(map_key, temp4));

                std::vector<int> temp5;
                temp5.push_back(iter_num);
                equilClasses.push_back(temp5);
                
            }

            goto LB1;


        }
        else
        {
            // multiple grps with max # of cases sharing same FGL - further selection based on length of sharing region
            for (std::map<int, std::vector<int>>::iterator iter1 = sharingGrpsRgV2.begin(); iter1 != sharingGrpsRgV2.end(); iter1++)
            {
                
                for (size_t i = 0; i < keys2keep.size(); i++)
                {
                    if (keys2keep[i] == iter1->first)
                    {
                        std::vector<int> temp1 = iter1->second;
                        shareLength.push_back(temp1[0]); // push length of contiguous FGL sharing
                    }

                }
                
            }


            if (!shareLength.empty())
            {
                maxLength = *max_element(shareLength.begin(), shareLength.end());
            }


            for (std::map<int, std::vector<int>>::iterator iter1 = sharingGrpsRgV2.begin(); iter1 != sharingGrpsRgV2.end(); iter1++)
            {
                for (size_t i = 0; i < keys2keep.size(); i++)
                {
                    if (keys2keep[i] == iter1->first)
                    {
                        std::vector<int> temp1 = iter1->second;
                        if (temp1[0] == maxLength)
                        {
                            keys2keep2.push_back(iter1->first);
                        }
                    }
                }
            }

            // record result using 0 el in keys2keep2

            std::string map_key;

            for (std::map<int, std::vector<int>>::iterator iter1 = sharingGrpsV2.begin(); iter1 != sharingGrpsV2.end(); iter1++)
            {

                int cur_key = iter1->first;

                if (keys2keep2[0] == cur_key)
                {

                    std::vector<int> temp1 = iter1->second;

                    for (size_t i = 3; i < temp1.size(); i++)
                    {
                        map_key += std::to_string(temp1[i]);
                    }
                }
            }

            if (sharingGrpsTotal.find(map_key) != sharingGrpsTotal.end()) // key exists
            {

                for (std::map<std::string, std::vector<int>>::iterator iter1 = sharingGrpsTotal.begin(); iter1 != sharingGrpsTotal.end(); iter1++)
                {
                    if (iter1->first == map_key)
                    {
                        std::vector<int> temp3 = iter1->second;

                        int el2edit = temp3[0];
                        equilClasses[el2edit].push_back(iter_num);

                    }
                }
            }
            else
            {
                mapKeyNum2++;

                std::vector<int> temp4;
                temp4.push_back(mapKeyNum2);

                for (std::map<int, std::vector<int>>::iterator iter1 = sharingGrpsV2.begin(); iter1 != sharingGrpsV2.end(); iter1++)
                {

                    if (iter1->first == keys2keep2[0])
                    {

                        std::vector<int> temp1 = iter1->second;
                    
                        for (size_t i = 3; i < temp1.size(); i++)
                        {
                            temp4.push_back(temp1[i]);
                        }
                    
                    }
                }

                sharingGrpsTotal.insert(make_pair(map_key, temp4));

                std::vector<int> temp5;
                temp5.push_back(iter_num);
                equilClasses.push_back(temp5);

            
            }

            goto LB1;


        }



    }
    else
    {

        if (sharingGrpsV2.size() == 1)
        {
            // record what you have
            std::string map_key;

            for (std::map<int, std::vector<int>>::iterator iter1 = sharingGrpsV2.begin(); iter1 != sharingGrpsV2.end(); iter1++)
            {

                std::vector<int> temp1 = iter1->second;

                for (size_t i = 3; i < temp1.size(); i++)
                {
                    map_key += std::to_string(temp1[i]);
                }
            }

            if (sharingGrpsTotal.find(map_key) != sharingGrpsTotal.end())
            {

                for (std::map<std::string, std::vector<int>>::iterator iter1 = sharingGrpsTotal.begin(); iter1 != sharingGrpsTotal.end(); iter1++)
                {
                    if (iter1->first == map_key)
                    {
                        std::vector<int> temp3 = iter1->second;

                        int el2edit = temp3[0];
                        equilClasses[el2edit].push_back(iter_num);

                    }
                }
            }
            else
            {
                mapKeyNum2++;

                std::vector<int> temp4;
                temp4.push_back(mapKeyNum2);

                for (std::map<int, std::vector<int>>::iterator iter1 = sharingGrpsV2.begin(); iter1 != sharingGrpsV2.end(); iter1++)
                {

                    std::vector<int> temp1 = iter1->second;
                   
                    for (size_t i = 3; i < temp1.size(); i++)
                    {
                        temp4.push_back(temp1[i]);
                    }
                    
                }

                sharingGrpsTotal.insert(make_pair(map_key, temp4));

                std::vector<int> temp5;
                temp5.push_back(iter_num);
                equilClasses.push_back(temp5);

                
               
            }
    
        }

    }

    LB1:;
}

/// identifies the core group of cases and iteration # to be used for phasing and haplotyping
std::vector<int> FGLshareGroup2::evalFGLshareLv3(std::map<std::string, std::vector<int>> &sharingGrpsTotal)
{

    std::vector<int> result1;// [0] = iteration # for phasing; [1]+ = IDs for cases from core group
    int max_el = -1;
    int iter_use;// iteration # to use for phasing and haplotyping

    std::vector<int> sizeEquilCl;

    for (size_t i = 0; i < equilClasses.size(); i++)
    {
        sizeEquilCl.push_back(equilClasses[i].size());
    }

    int maxNumIter = -1;
    if (!sizeEquilCl.empty())
    {
        maxNumIter = *max_element(sizeEquilCl.begin(), sizeEquilCl.end());
    }

    for (size_t i = 0; i < equilClasses.size(); i++)
    {
        if (equilClasses[i].size() == static_cast<size_t>(maxNumIter))
        {
            max_el = i;
            iter_use = equilClasses[i][0];
            goto LB1;
        }
    }
    LB1:;


    for (auto iter1 = sharingGrpsTotal.begin(); iter1 != sharingGrpsTotal.end(); iter1++)
    {
        vector<int> inVect = (*iter1).second;

        if (inVect[0] == max_el)
        {
            result1.push_back(iter_use);
           
            for (size_t j = 1; j < inVect.size(); j++)
            {
                result1.push_back(inVect[j]);
            }
           
            goto LB2;
        }
    }
    LB2:;

    return result1;
}


class MatrixHaplo // subject chr_n, and pos_n start with 0
{
public:

    MatrixHaplo(int dim3, std::vector<int> *current_wind_bps, Pedigree_info *current_ped_ptr, int chrNumIn);

    ~MatrixHaplo();
    
    void setVal(int, int, int);

    void setHeader();

    void printMatrixLine(int);

    void printOutResult(ofstream *haploStream2);

    void printOutResultVCF(ofstream *haploStream3);


private:
    int num_pos;
    std::vector<int> *cur_window_bps = nullptr;
    Pedigree_info *ped1 = nullptr;
    int num_chrs;
    int **matrixHaplo = nullptr;
    int chrNum = 0;
    std::vector<std::vector<std::string>> storage;
};

MatrixHaplo::MatrixHaplo(int dim3, std::vector<int> *current_wind_bps,
Pedigree_info *current_ped_ptr, int chrNumIn) : num_pos{dim3},
cur_window_bps{current_wind_bps}, ped1{current_ped_ptr}, chrNum{chrNumIn}
{

    num_chrs = (ped1->getNumWGS()) * 2;

    matrixHaplo = new int *[num_chrs];
    for (int i = 0; i < num_chrs; i++)
    {

        matrixHaplo[i] = new int[num_pos];
        for (int k = 0; k < num_pos; k++)
        {
            matrixHaplo[i][k] = -1;
        }
    }

    storage.resize(num_chrs);
}

MatrixHaplo::~MatrixHaplo()
{

    for (int i = 0; i < num_chrs; i++)
    {
        delete[] matrixHaplo[i];
    }
    delete[] matrixHaplo;
}

void MatrixHaplo::setVal(int chr_n_3D, int pos_n, int al) // all variables start with 0
{
    matrixHaplo[chr_n_3D][pos_n] = al;
}

void MatrixHaplo::setHeader()
{
    std::string chr = std::to_string(chrNum);
    std::string first = std::to_string(cur_window_bps->front());
    std::string last = std::to_string(cur_window_bps->back());
    std::string window_range = chr + ":" + first + "-" + last;

    for (size_t i = 0; i < ped1->cur_ped.size(); i++)
    {

        if (ped1->cur_ped[i].wgs_data == 1)
        {
            storage[ped1->cur_ped[i].mat_3D_pos].push_back(window_range);
            storage[ped1->cur_ped[i].pat_3D_pos].push_back(window_range);

            storage[ped1->cur_ped[i].mat_3D_pos].push_back(ped1->cur_ped[i].subject_id_orig + "_0");
            storage[ped1->cur_ped[i].pat_3D_pos].push_back(ped1->cur_ped[i].subject_id_orig + "_1");

            printMatrixLine(ped1->cur_ped[i].mat_3D_pos);
            printMatrixLine(ped1->cur_ped[i].pat_3D_pos);
        }
    }
}

void MatrixHaplo::printMatrixLine(int chr_num_3D) // all variables start with 0
{

    int k = chr_num_3D;

    std::string temp1;
    for (int j = 0; j < num_pos; j++)
    {
        temp1.append(to_string(matrixHaplo[k][j]));
    }
    storage[k].push_back(temp1); // k should be mat_3D_pos or pat_3D_pos
}

void MatrixHaplo::printOutResult(ofstream *haploStream2)
{

    setHeader();

    for (size_t i = 0; i < storage.size(); i++)
    {
        for (size_t j = 0; j < storage[i].size(); j++)
        {

            *haploStream2 << storage[i][j] << " ";
        }

        *haploStream2 << std::endl;
    }
}


class Phase2
{

public:
    Phase2(FGL_matrix *, DenseMarkerGenos *, Pedigree_info *, MatrixHaplo *, int, int, int, ofstream *);

    void initializePhasedGeno();

    void phaseGenos();

    int assignGenoToFGL();

    int assignHomGenoToFGL();

    int checkFGLGenoConsistency();

    int resolveHets();

    int countResolvedFGL();

    void assignHomGenos2Chrs();
    
    void sendPhasedGenoToHaploMatrix();

    void sendPhasedGenoToHaploMatrix2();

private:
    FGL_matrix *matrix_FGL = nullptr;

    DenseMarkerGenos *dense_genotypes = nullptr;

    Pedigree_info *ped1 = nullptr;

    MatrixHaplo *haplo_window = nullptr;

    ofstream *consistencyStream = nullptr;

    int var_bp;

    int var_pos;

    int cur_iter_num = -1;

    std::vector<std::vector<std::string>> phased_geno;
};

Phase2::Phase2(FGL_matrix *cur_matrix_FGL, DenseMarkerGenos *cur_dense_genotypes,
Pedigree_info *cur_ped1, MatrixHaplo *haplo_window_ptr, int cur_var_bp, int cur_var_pos,
int iter_num, ofstream *consistencyStreamIn) : matrix_FGL{cur_matrix_FGL}, dense_genotypes{cur_dense_genotypes},
ped1{cur_ped1}, haplo_window{haplo_window_ptr}, consistencyStream{consistencyStreamIn}, var_bp{cur_var_bp}, var_pos{cur_var_pos},
cur_iter_num{iter_num}
{
    this->initializePhasedGeno();
}

/// initializes data structure to store FGL assigned genotype alleles
void Phase2::initializePhasedGeno()
{
    int num_chrs = (ped1->getNumFounders()) * 2; // number of FGLs
    phased_geno.resize(num_chrs);
}

/// main function for this class
void Phase2::phaseGenos() 
{

    int outcome1 = -1;

    outcome1 = this->assignGenoToFGL();

    if (outcome1 == 1 || outcome1 == 2)
    {
        this->sendPhasedGenoToHaploMatrix();
    }

    else if (outcome1 == 0)
    {
        *consistencyStream << var_bp << std::endl;
    }

}

/// sub-function for this class, performs phasing
int Phase2::assignGenoToFGL()
{

    int consistency1 = -1;
    int check1 = -1;
    int resolved_fgl = 0;

    check1 = this->assignHomGenoToFGL();

    if (check1 == 0)
    {
        return 2; // exit, cannot phase because no homozygous genotypes
    }

    consistency1 = this->checkFGLGenoConsistency(); // if consistency == 1 then no problems

    if (consistency1 == 1)
    {
    LB1:

        int consistency2 = -1;
        consistency2 = this->resolveHets();

        if (consistency2 == 1)
        {
            int consistency3 = -1;
            consistency3 = this->checkFGLGenoConsistency();

            if (consistency3 == 1)
            {
                int temp_res1 = this->countResolvedFGL();

                if (temp_res1 > resolved_fgl) // keep resolving fgl until reach plateau
                {
                    resolved_fgl = temp_res1;

                    goto LB1;
                }
                else
                {
                    return 1; // successful end
                }
            }
            else
            {
                this->assignHomGenos2Chrs();
                this->sendPhasedGenoToHaploMatrix2();
                return 0;
            }
        }
        else
        {
            this->assignHomGenos2Chrs();
            this->sendPhasedGenoToHaploMatrix2();
            return 0;
        }
    }
    else
    {
        this->assignHomGenos2Chrs();
        this->sendPhasedGenoToHaploMatrix2();
        return 0;
    }
}

/// checks if there are any homozygous genotypes and assigns their alleles to FGLs
int Phase2::assignHomGenoToFGL()
{

    int check1 = 0; // to check if we have at least one homozygous genotype

    for (size_t i = 0; i < ped1->cur_ped.size(); i++)
    {

        if (ped1->cur_ped[i].wgs_data == 1)
        {
            std::vector<std::string> result1;

            result1 = dense_genotypes->getGeno(var_bp, i);

            if ((result1[0].compare(result1[1]) == 0) &&
            ((result1[0].compare("0") != 0) && (result1[1].compare("0") != 0)))
            {
                check1++;
                std::vector<int> result2;
                result2 = matrix_FGL->getFGL(cur_iter_num, i);

                phased_geno[result2[0] - 1].push_back(result1[0]);
                phased_geno[result2[1] - 1].push_back(result1[1]);
            }
        }
    }

    if (check1 > 0)
    {
        return 1; 
    }
    else
    {
        return 0; // if no homozygous genotypes are present
    }
}

/// checks for consistency assignment of alleles from homozygous genotypes to FGLs
int Phase2::checkFGLGenoConsistency()
{

    for (size_t it1 = 0; it1 < phased_geno.size(); it1++)
    {

        if (phased_geno[it1].size() > 1)
        {
            std::string al1 = phased_geno[it1][0];

            for (size_t it2 = 1; it2 < phased_geno[it1].size(); it2++)
            {

                if (phased_geno[it1][it2].compare(al1) == 0)
                {
                    //  do nothing
                }
                else
                {
                    return 0; // Problem, not consistent
                }
            }
        }
    }

    return 1; // all consistent
}

/// assigns alleles of heterozygous genotypes to FGLs 
int Phase2::resolveHets()
{

    for (size_t i = 0; i < ped1->cur_ped.size(); i++)
    {

        if (ped1->cur_ped[i].wgs_data == 1)
        {
            std::vector<std::string> result1;

            result1 = dense_genotypes->getGeno(var_bp, i);

            if ((result1[0].compare(result1[1]) != 0) &&
            ((result1[0].compare("0") != 0) && (result1[1].compare("0") != 0)))
            {

                std::vector<int> result2;
                result2 = matrix_FGL->getFGL(cur_iter_num, i);

                if (!phased_geno[result2[0] - 1].empty() &&
                !phased_geno[result2[1] - 1].empty()) // just check for consistency
                {

                    if (phased_geno[result2[0] - 1][0].compare(phased_geno[result2[1] - 1][0]) == 0)
                    {
                        return 0; // Problem, inconsistent
                    }
                    
                }

                else if (!phased_geno[result2[0] - 1].empty() && phased_geno[result2[1] - 1].empty())
                {
                   
                    if (phased_geno[result2[0] - 1][0].compare("1") == 0)
                    {
                        phased_geno[result2[1] - 1].push_back("2");
                    }
                    else
                    {
                        phased_geno[result2[1] - 1].push_back("1");
                    }
                }

                else if (!phased_geno[result2[1] - 1].empty() && phased_geno[result2[0] - 1].empty())
                {
                    
                    if (phased_geno[result2[1] - 1][0].compare("1") == 0)
                    {
                        phased_geno[result2[0] - 1].push_back("2");
                    }
                    else
                    {
                        phased_geno[result2[0] - 1].push_back("1");
                    }
                }
            }
        }
    }

    return 1; // still need to check for consistency later
}

/// counts the number of FGLs which got allele assignment after phasing
int Phase2::countResolvedFGL()
{

    int resolved_fgl = 0;

    for (size_t it1 = 0; it1 < phased_geno.size(); it1++)
    {

        if (phased_geno[it1].size() > 0)
        {
            resolved_fgl++;
        }
    }

    return resolved_fgl;
}

/// assigns alleles of homozygous genotypes to FGLs if phasing is not possible due to inconsistencies
void Phase2::assignHomGenos2Chrs()
{

    phased_geno.clear(); 
    this->initializePhasedGeno(); 

    for (size_t i = 0; i < ped1->cur_ped.size(); i++)
    {

        if (ped1->cur_ped[i].wgs_data == 1)
        {
            std::vector<std::string> result1;

            result1 = dense_genotypes->getGeno(var_bp, i);

            if ((result1[0].compare(result1[1]) == 0) &&
            ((result1[0].compare("0") != 0) && (result1[1].compare("0") != 0)))
            {

                std::vector<int> result2;
                result2 = matrix_FGL->getFGL(cur_iter_num, i);

                phased_geno[result2[0] - 1].push_back(result1[0]);
                phased_geno[result2[1] - 1].push_back(result1[1]);
                
            }

        }

    }

    this->sendPhasedGenoToHaploMatrix2();

}

/// stores phased genotype data
void Phase2::sendPhasedGenoToHaploMatrix()
{ 

    for (size_t i = 0; i < ped1->cur_ped.size(); i++)
    {

        if (ped1->cur_ped[i].wgs_data == 1)
        {
            std::vector<std::string> result1;
            result1 = dense_genotypes->getGeno(var_bp, i);
            std::vector<int> result2;
            result2 = matrix_FGL->getFGL(cur_iter_num, i);

            if ((result1[0].compare("0") != 0) && (result1[1].compare("0") != 0))
            {

                if (!phased_geno[result2[0] - 1].empty())
                {

                    haplo_window->setVal(ped1->cur_ped[i].mat_3D_pos,
                    var_pos, stoi(phased_geno[result2[0] - 1][0]));
                }

                else

                {

                    haplo_window->setVal(ped1->cur_ped[i].mat_3D_pos,
                    var_pos, 3);
                }

                if (!phased_geno[result2[1] - 1].empty())
                {

                    haplo_window->setVal(ped1->cur_ped[i].pat_3D_pos,
                    var_pos, stoi(phased_geno[result2[1] - 1][0]));
                }

                else

                {

                    haplo_window->setVal(ped1->cur_ped[i].pat_3D_pos,
                    var_pos, 3);
                }
            }

            else

            {

                if (!phased_geno[result2[0] - 1].empty())
                {

                    haplo_window->setVal(ped1->cur_ped[i].mat_3D_pos,
                    var_pos, stoi(phased_geno[result2[0] - 1][0]));
                }

                else

                {

                    haplo_window->setVal(ped1->cur_ped[i].mat_3D_pos,
                    var_pos, 0);
                }

                if (!phased_geno[result2[1] - 1].empty())
                {

                    haplo_window->setVal(ped1->cur_ped[i].pat_3D_pos,
                    var_pos, stoi(phased_geno[result2[1] - 1][0]));
                }

                else

                {

                    haplo_window->setVal(ped1->cur_ped[i].pat_3D_pos,
                    var_pos, 0);
                }
            }
        }
    }
}

/// stores phased genotype data (only alleles from homozygous genotype) after inconsistency was detected
void Phase2::sendPhasedGenoToHaploMatrix2()
{ 

    for (size_t i = 0; i < ped1->cur_ped.size(); i++)
    {

        if (ped1->cur_ped[i].wgs_data == 1)
        {
            std::vector<std::string> result1;
            result1 = dense_genotypes->getGeno(var_bp, i);
            std::vector<int> result2;
            result2 = matrix_FGL->getFGL(cur_iter_num, i);

            if (((result1[0].compare("0") != 0) && (result1[1].compare("0") != 0)) &&
            ((result1[0].compare(result1[1]) == 0)))
            {

                haplo_window->setVal(ped1->cur_ped[i].mat_3D_pos,
                var_pos, stoi(result1[0]));
                
                haplo_window->setVal(ped1->cur_ped[i].pat_3D_pos,
                var_pos, stoi(result1[1]));
                
            }

            else if (((result1[0].compare("0") != 0) && (result1[1].compare("0") != 0)) &&
            ((result1[0].compare(result1[1]) != 0)))
            {

                haplo_window->setVal(ped1->cur_ped[i].mat_3D_pos, var_pos, 7);

                haplo_window->setVal(ped1->cur_ped[i].pat_3D_pos, var_pos, 7);


            }

            else if ((result1[0].compare("0") == 0) && (result1[1].compare("0") == 0))
            {

                haplo_window->setVal(ped1->cur_ped[i].mat_3D_pos, var_pos, 0);

                haplo_window->setVal(ped1->cur_ped[i].pat_3D_pos, var_pos, 0);

            }

           
            
        }
    }
}


class IntegrationHpSh
{
public:

    IntegrationHpSh(int, int, Pedigree_info *, ofstream *, ofstream *, Parameters *, DenseMarkerPos *);
    
    ~IntegrationHpSh();

    void setValue(int, std::vector<std::vector<int>> &);

    void setValue2(int, std::vector<int> &);

    void setValue3(int, std::vector<int> &);

    std::vector<int> readHaploShareData(int, string);

    void record2ndSharedHaplo(std::string, std::string, std::vector<std::string> &);

    int integrateHpSh();

    std::vector<int> resolveDoubleHaploSharing();

    void updateShareHaploSeq(std::vector<int> &, std::string);

    void printFileHeader();

    void printHaploSharing();

private:
    int **matrixHpShareAll = nullptr;
    int num_chrs;
    int num_wind;
    Pedigree_info *ped1 = nullptr;
    Parameters *pars = nullptr;
    DenseMarkerPos *dense_m_pos = nullptr;
    ofstream *haploSharingStream = nullptr;
    ofstream *sharedHaploStreamFinal = nullptr;
    std::vector<std::vector<std::vector<string>>> second_shared_haplo;
    std::vector<std::vector<std::vector<int>>> second_haplo_sharing; // V1 - case #; V2[0] = window#; V2[1] = 3D_chr for subjects from the candidate group who share a haplotype (unambig); V2[2] = 3D_chr for subjects from the candidate group who share a haplotype (ambig); v2[3] = 3D_chr other subjects with a copy of risk haplotype
    std::vector<std::string> candidate_grp;
};

IntegrationHpSh::IntegrationHpSh(int numSubWGS, int n_wind, Pedigree_info *current_ped_ptr,
ofstream *haploSharingStreamIn, ofstream *sharedHaploStreamFinalIn, Parameters *current_pars,
DenseMarkerPos *dense_marker_pos) : matrixHpShareAll{nullptr}, num_chrs{0}, num_wind{n_wind}, 
ped1{current_ped_ptr}, pars{current_pars}, dense_m_pos{dense_marker_pos}, haploSharingStream{haploSharingStreamIn}, 
sharedHaploStreamFinal{sharedHaploStreamFinalIn}
{

    // this matrix shows whether a chr is involved in haplotype sharing in a genomic window
    num_chrs = numSubWGS * 2;
    matrixHpShareAll = new int *[num_wind];
    for (int i = 0; i < num_wind; i++)
    {

        matrixHpShareAll[i] = new int[num_chrs];
        for (int k = 0; k < num_chrs; k++)
        {
            matrixHpShareAll[i][k] = -1;
        }
    }
}

IntegrationHpSh::~IntegrationHpSh()
{

    for (int i = 0; i < num_wind; i++)
    {
        delete[] matrixHpShareAll[i];
    }
    delete[] matrixHpShareAll;
}

/// records results of haplotype sharing determination for a genomic window
void IntegrationHpSh::setValue(int cur_wind, std::vector<std::vector<int>> &shared_chr_3D)
{

    /**k values are 0 to 3, which will translate into following values for:
     * sharing grp 1: 
     * 1 = shared chr, core set; 
     * 2 = ambig shared chr, core set;
     * 3 = shared chr, others (not from core set);
     * 4 = ambig shared chr, others;
     * sharing grp 2:
     * 5 = shared chr, core set; 
     * 6 = ambig shared chr with double sharing, core set:
     * 7 = shared chr, others (not from core set);
     * 8 = ambig shared chr with double sharing, others;
     */
    for (size_t k = 0; k < shared_chr_3D.size(); k++) 
    {
        for (size_t j = 0; j < shared_chr_3D[k].size(); j++)
        {

            matrixHpShareAll[cur_wind - 1][shared_chr_3D[k][j]] = k + 1;
        }
    }
}

/// records results of haplotype sharing determination for a genomic window
void IntegrationHpSh::setValue2(int cur_wind, std::vector<int> &candidate_grp3Dchrs)
{

    for (size_t k = 0; k < candidate_grp3Dchrs.size(); k++) 
    {

        if (matrixHpShareAll[cur_wind - 1][candidate_grp3Dchrs[k]] == -1)
        {
            matrixHpShareAll[cur_wind - 1][candidate_grp3Dchrs[k]] = 5;
        }

        else if (matrixHpShareAll[cur_wind - 1][candidate_grp3Dchrs[k]] == 2)
        {
            matrixHpShareAll[cur_wind - 1][candidate_grp3Dchrs[k]] = 6;
        }

    }
}

/// records results of haplotype sharing determination for a genomic window
void IntegrationHpSh::setValue3(int cur_wind, std::vector<int> &others3Dchrs)
{

    for (size_t k = 0; k < others3Dchrs.size(); k++) 
    {

        if (matrixHpShareAll[cur_wind - 1][others3Dchrs[k]] == -1)
        {
            matrixHpShareAll[cur_wind - 1][others3Dchrs[k]] = 7;
        }

        else if (matrixHpShareAll[cur_wind - 1][others3Dchrs[k]] == 4)
        {
            matrixHpShareAll[cur_wind - 1][others3Dchrs[k]] = 8;
        }

    }
}

/// returns haplotype sharing status for both chrs of a subject in a genomic window (could be -1 or any whole number between 1 and 8)
std::vector<int> IntegrationHpSh::readHaploShareData(int cur_wind, string subjectID)
{

    std::vector<int> result;

    int subjectIndex = ped1->getSubjectIndex(subjectID);

    result.push_back(matrixHpShareAll[cur_wind - 1][ped1->cur_ped[subjectIndex].mat_3D_pos]);
    result.push_back(matrixHpShareAll[cur_wind - 1][ped1->cur_ped[subjectIndex].pat_3D_pos]);

    return result;
}

/// records 2nd shared haplotype sequence
void IntegrationHpSh::record2ndSharedHaplo(std::string cur_wind_el_str, std::string windBpBounds, std::vector<std::string> &secondSharedHaploSeq)
{

    std::vector<std::vector<string>> temp1(3);
    std::vector<string> temp2;
    std::vector<string> temp3;

    temp2.push_back(cur_wind_el_str); // this is window# - 1
    temp3.push_back(windBpBounds);

    temp1[0].insert(temp1[0].end(), temp2.begin(), temp2.end());
    temp1[1].insert(temp1[1].end(), temp3.begin(), temp3.end());
    temp1[2].insert(temp1[2].end(), secondSharedHaploSeq.begin(), secondSharedHaploSeq.end());

    second_shared_haplo.push_back(temp1);
}

/// main function for the class that integrates entire workflow
int IntegrationHpSh::integrateHpSh()
{

    int result = -1;

    std::vector<int> change2secondHpShare = this->resolveDoubleHaploSharing();// contains window # where change has to be made

    if (!change2secondHpShare.empty())
    {
        result = 1;

        this->updateShareHaploSeq(change2secondHpShare, pars->output_dir_path);

    }

    this->printFileHeader();

    this->printHaploSharing();

    return result;
}

/// determines which shared haplotype is the risk haplotype
std::vector<int> IntegrationHpSh::resolveDoubleHaploSharing()
{

    std::vector<int> result1; // push window numbers where 2nd shared haplotype sequence should be printed into risk haplotype sequence file
    std::vector<int> wind2examine;// 1st and last genomic window # in ROI
    std::vector<double> linkage_boundaries;// ROI boundaries in cM
    std::vector<std::vector<std::string>> subjectsWithSharedHaplo(2); // 0 = ids with maternal chr; 1 = ids with parternal chr
 
    linkage_boundaries = pars->linkage_region;

    wind2examine = dense_m_pos->getROIboundWindNum(linkage_boundaries);


    for (size_t i = 0; i < candidate_grp.size(); i++)
    {
       
        int m_chr_count = 0;
        int p_chr_count = 0;
        

        for (int k = wind2examine[0]; k <= wind2examine[1]; k++) 
        {
            std::vector<int> subjectHaploshare =  this->readHaploShareData(k, candidate_grp[i]);


            if (subjectHaploshare[0] == 1 || subjectHaploshare[0] == 5)
            {
                m_chr_count++;
            }

            if (subjectHaploshare[1] == 1 || subjectHaploshare[1] == 5)
            {
                p_chr_count++;
            }

        }


        if (m_chr_count >= p_chr_count)
        {

            subjectsWithSharedHaplo[0].push_back (candidate_grp[i]);
            
        }
        else
        {

            subjectsWithSharedHaplo[1].push_back (candidate_grp[i]);
            
        }


    }

// now determine which shared haplotype sequence to use for risk haplotype

    for (int k = wind2examine[0]; k <= wind2examine[1]; k++) 
    {

        int count_1_2 = 0;
        int count_5_6 = 0;


        for (size_t i = 0; i < subjectsWithSharedHaplo.size(); i++)
        {

            for (size_t j = 0; j < subjectsWithSharedHaplo[i].size(); j++)
            {

                std::vector<int> subjectHaploshare =  this->readHaploShareData(k, subjectsWithSharedHaplo[i][j]);

                if (i == 0) // for maternal cases
                {
                    if (subjectHaploshare[0] == 1 || subjectHaploshare[0] == 2)
                    {
                        count_1_2++;
                    }
                    else if (subjectHaploshare[0] == 5 || subjectHaploshare[0] == 6)
                    {
                        count_5_6++;
                    }

                }
                else // for parternal cases
                {
                    if (subjectHaploshare[1] == 1 || subjectHaploshare[1] == 2)
                    {
                        count_1_2++;
                    }
                    else if (subjectHaploshare[1] == 5 || subjectHaploshare[1] == 6)
                    {
                        count_5_6++;
                    }


                }

            }

        }

        if (count_5_6 > count_1_2)
        {
            result1.push_back (k);
        }
       

    }


    return result1;
}

/// outputs risk haplotype sequence into a file
void IntegrationHpSh::updateShareHaploSeq(std::vector<int> &change2secondHpShare, std::string outputDir)
{

    std::map<std::string, int> windowIndex; // window# - 1

    for (size_t k = 0; k < change2secondHpShare.size(); k++)
    {

        for (size_t j = 0; j < second_shared_haplo.size(); j++)
        {

            if (change2secondHpShare[k] - 1 == std::stoi(second_shared_haplo[j][0][0])) // window# - 1
            {
                windowIndex.insert(std::make_pair(second_shared_haplo[j][1][0], j));// bp position boundaries of a genomic window and el # of second_shared_haplo with haplo sequence for that window 
            }
        }
    }

   
    std::string sharedHaplotypeTempFileName = outputDir + "/shared_haplotype_temp.txt";
    std::ifstream inSharedHaplotypeTempFile;
    inSharedHaplotypeTempFile.open(sharedHaplotypeTempFileName, std::ios::in);

    if (!inSharedHaplotypeTempFile)
    {
        cerr << "shared_haplotype_temp.txt file cannot be opened" << std::endl;
        exit(1);
    }

    std::string line;

    if (inSharedHaplotypeTempFile.is_open())
    {
        while (getline(inSharedHaplotypeTempFile, line))
        {

            std::string windBpBound;
            std::stringstream ss(line);
            ss >> windBpBound;

            std::map<std::string, int>::iterator it;
            it = windowIndex.find(windBpBound);


            if (it != windowIndex.end())
            {

                int index = windowIndex.find(windBpBound)->second;

                std::string tempHaploSeq;

                for (size_t g = 0; g < second_shared_haplo[index][2].size(); g++)
                {
                    tempHaploSeq.append(second_shared_haplo[index][2][g]);
                }

                *sharedHaploStreamFinal << second_shared_haplo[index][1][0] << " " << tempHaploSeq << std::endl;
            }

            else
            {
                *sharedHaploStreamFinal << line << std::endl;
            }
        }
    }

    inSharedHaplotypeTempFile.close();
}


/// prints out a header for haplotype sharing pattern file
void IntegrationHpSh::printFileHeader()
{

    *haploSharingStream << "Genomic_window";

    for (int k = 0; k < num_chrs; k++)
    {

        for (size_t i = 0; i < ped1->cur_ped.size(); i++)
        {

            if (ped1->cur_ped[i].mat_3D_pos == k ||
                ped1->cur_ped[i].pat_3D_pos == k)
            {
                *haploSharingStream << " " << ped1->cur_ped[i].subject_id_orig;
            }
        }
    }

    *haploSharingStream << std::endl;
}


/// prints haplotype sharing results to a file; for each chr of subjects with WGS data: 0 = no sharing; 1 = sharing between cases from the candidate group; 2 = sharing of others with cases from the candidate group
void IntegrationHpSh::printHaploSharing()
{

    for (int i = 0; i < num_wind; i++)
    {

        int wind_num = i + 1;
        *haploSharingStream << wind_num;

        for (int k = 0; k < num_chrs; k++)
        {

            if (matrixHpShareAll[i][k] == -1)
            {
                *haploSharingStream << " " << 0;
            }

            else if (matrixHpShareAll[i][k] == 1 || matrixHpShareAll[i][k] == 2 || \
            matrixHpShareAll[i][k] == 5 || matrixHpShareAll[i][k] == 6)
            {
                *haploSharingStream << " " << 1;
            }

            else if (matrixHpShareAll[i][k] == 3 || matrixHpShareAll[i][k] == 4 || \
            matrixHpShareAll[i][k] == 7 || matrixHpShareAll[i][k] == 8)
            {
                *haploSharingStream << " " << 2;
            }

        }

        *haploSharingStream << std::endl;
    }
}


class HaploShare
{
public:
    HaploShare(int, int, Pedigree_info *, IntegrationHpSh *, Parameters &,
    DenseMarkerGenos *, std::vector<int> *, ofstream *);

    ~HaploShare();
    
    void setHaploData(int, const std::vector<int> &);

    void printHaploData(int, int);

    void makeAOA1();

    void makeAOA2();

    void makeAOA3();

    int selectRoute();

    std::vector<int> getColNumComb(int);

    std::vector<std::vector<std::string>> convert3DChr2Sub(std::vector<int> &);

    void search4HaploShare(int);

    void testColSubstAOA(int, std::vector<int> &);

    void analyzeHpShP1();

    void initAOA4Subst(int);

    void checkCandidateGrpMatch(std::vector<int> &, std::vector<int> &);

    int check2ndSharedHaploPossibility();

    void determine2ndSharedHaploSeq();

    std::vector<int> determineOthersWhoShare2ndHaplo(std::vector<int> &);

    void setCurSharingInfo(std::vector<int> &, std::vector<int> &, std::vector<int> &, std::vector<int> &);

    void selectBestCurHaploShare();

    std::vector<int> selectCurHaploSharesStep1();

    std::vector<int> selectCurHaploSharesStep2(std::vector<int> &);

    std::vector<int> selectCurHaploSharesStep3(std::vector<int> &);

    void selectCurHaploSharesStep4(std::vector<int> &);

    std::vector<int> checkHpSharingMatrixTemp1();

    std::vector<int> getVarPos2VarifyAOA1(std::vector<int> &);

    std::vector<std::vector<int>> checkGeno4ProposedAl(std::vector<int> &);

    std::vector<std::vector<int>> checkGeno4ProposedAlOthers(std::vector<std::string> &, int);

    int getNumChrsShared(std::string &);

    std::vector<std::string> getSubsWithBothChrs(std::vector<string> &);

    std::vector<int> get3Dchrs4candidate_grp();

    std::vector<int> detectAmbigSharedHaplosPresence();

    void evalAmbigSharedHaplosPresence(std::vector<int> &);

private:
    int **matrixHaploData = nullptr;
    int num_pos;
    int num_chrs;
    int cur_wind; // starts from 1
    Pedigree_info *ped1 = nullptr;
    IntegrationHpSh *global_hp_share = nullptr;
    DenseMarkerGenos *dense_genotypes = nullptr;
    std::vector<int> *cur_window_bps = nullptr;

    int num_rows_aoa = 0;
    int chr = 0;
    int **matrixAOA1 = nullptr;
    int **matrixAOA2 = nullptr;
    int **matrixAOA3 = nullptr;
    int **matrixTemp1 = nullptr;

    std::vector<std::vector<int>> aoa_chrs;          // element of 1st vector indicates haplo group number and elements of 2nd vector contain chr origin of subjects contributing to haplo group
    std::vector<std::vector<std::string>> aoa_ids;   // element of 1st vector indicates haplo group number and elements of 2nd vector contain ids of subjects contributing to haplo group
    std::vector<std::vector<std::string>> aoa_2chrs; // element of 1st vector indicates haplo group number and elements of 2nd vector contain ids of subjects contributing to haplo group both of their chrs

    std::vector<std::vector<std::vector<int>>> cur_sharing; // V1 = case count #; V2[0] = substitution combination or -1; V2[1] = row #s for sharing group; V2[2] = suggested alleles for unphased variants or -1; V2[3] = subject index of those from core group having only one chr IBS with shared haplotype, V2[4] = subject index of those from non candidate group possessing at least a copy of shared haplotype, V2[5] = subject index of those from non candidate group having only one chr IBS with shared haplotype

    std::vector<std::string> candidate_grp;

    ofstream *sharedHaploStream = nullptr;

    std::vector<int> sharedHaploOut;
};

HaploShare::HaploShare(int window_size, int cur_wind_in, Pedigree_info *current_ped_ptr,
IntegrationHpSh *global_hp_share_ptr, Parameters &current_pars, DenseMarkerGenos *cur_dense_genotypes,
std::vector<int> *current_wind_bps, ofstream *sharedHaploStreamIn) : matrixHaploData{nullptr},                     
num_pos{window_size}, num_chrs{0}, cur_wind{cur_wind_in}, ped1{current_ped_ptr},                      
global_hp_share{global_hp_share_ptr}, dense_genotypes{cur_dense_genotypes}, cur_window_bps{current_wind_bps},            
candidate_grp{current_pars.candidate_group}, sharedHaploStream{sharedHaploStreamIn}    
{

    chr = current_pars.chr;

    num_chrs = (ped1->getNumWGS()) * 2;

    matrixHaploData = new int *[num_chrs];
    for (int i = 0; i < num_chrs; i++)
    {

        matrixHaploData[i] = new int[num_pos];
        for (int k = 0; k < num_pos; k++)
        {
            matrixHaploData[i][k] = -1;
        }
    }
}

HaploShare::~HaploShare()
{

    for (int i = 0; i < num_chrs; i++)
    {
        delete[] matrixHaploData[i];
    }
    delete[] matrixHaploData;

    for (int i = 0; i < num_rows_aoa; i++)
    {
        delete[] matrixAOA1[i];
    }
    delete[] matrixAOA1;

    for (int i = 0; i < num_rows_aoa; i++)
    {
        delete[] matrixAOA2[i];
    }
    delete[] matrixAOA2;

    for (int i = 0; i < num_rows_aoa; i++)
    {
        delete[] matrixAOA3[i];
    }
    delete[] matrixAOA3;

    for (int i = 0; i < num_rows_aoa; i++)
    {
        delete[] matrixTemp1[i];
    }
    delete[] matrixTemp1;
}

/// reading haplotype sequences into data structure
void HaploShare::setHaploData(int chr_num_3D, const std::vector<int> &haplo_seq)
{

    for (size_t k = 0; k < haplo_seq.size(); k++)
    {
        matrixHaploData[chr_num_3D][k] = haplo_seq[k];
    }
}

/// printing haplotype sequences stored in HaploData structure
void HaploShare::printHaploData(int chr_num, int var_num)
{
    
    for (int k = 0; k < chr_num; k++)
    {
        for (int j = 0; j < var_num; j++)
        {
            std::cout << matrixHaploData[k][j] << "  ";
        }
        std::cout << std::endl;
    }

}


/// this is the main function for this class which executes other class functions to determine haplotype sharing
void HaploShare::analyzeHpShP1()
{

    this->makeAOA1();

    this->makeAOA2();

    this->makeAOA3();

    int route_num = this->selectRoute(); // outcomes are 0, 1, or 2

    this->search4HaploShare(route_num);

    this->selectBestCurHaploShare();

    if (!cur_sharing.empty())
    {

        int secondSharingRes = this->check2ndSharedHaploPossibility(); // if -1 then there is no 2nd sharing; if 1 then 2nd sharing exists


        if (secondSharingRes == 1)
        {
            this->determine2ndSharedHaploSeq();
        }
        else
        {
            
            // check if 2nd ambiguouse sharing exist
            std::vector<int> result1 = this->detectAmbigSharedHaplosPresence();

            if(result1[0] != -1)
            {
                // compare sequences of ambig shared haplotypes to check if they are the same
                this->evalAmbigSharedHaplosPresence(result1);
            }

        }
    }

}

/// creating AOA1 matrix of unique haplotype sequences, where REF and ALT alleles are coded as 1 and 2, and missing or unphasable are coded as 0
void HaploShare::makeAOA1()
{

    std::vector<std::vector<int>> arrayHaplo_count;

    NumberList listIter;
    for (int b = 0; b < num_chrs; b++)
    {
        listIter.appendNode(b);
    }

    while (listIter.countListNodes() >= 2)
    {

        std::vector<int> values_delete;
        std::vector<int> temp_same_haplo_chr;
        int st = listIter.ptrNodeElement(0)->value;
        temp_same_haplo_chr.push_back(st);

        for (int it = 1; it < listIter.countListNodes(); it++) // element number in the list
        {

            int i = listIter.ptrNodeElement(it)->value;

            int count2 = 0; // iterates over variant position number within haplo

            for (int j = 0; j < num_pos; j++) // variant position #
            {
                if (matrixHaploData[st][j] != matrixHaploData[i][j])
                {

                    goto LABEL1;
                }
                else
                {
                    count2++;
                }
            }

            if (count2 == num_pos)
            {

                temp_same_haplo_chr.push_back(i);

                values_delete.push_back(i);
            }

        LABEL1:; // go to next chr
        }

        if (!values_delete.empty())
        {

            for (std::vector<int>::iterator iter1 = values_delete.begin();
                 iter1 != values_delete.end(); ++iter1)
            {
                listIter.deleteNode(*iter1);
            }
        }

        listIter.deleteNode(st);

        arrayHaplo_count.push_back(temp_same_haplo_chr);

    } // end of while statement


    if(listIter.countListNodes())
    {
        std::vector<int> temp_same_haplo_chr;
        int i = listIter.ptrNodeElement(0)->value;
        temp_same_haplo_chr.push_back(i);
        arrayHaplo_count.push_back(temp_same_haplo_chr);
    }


    // sorting to put haplotype sequence present on the largest # of chrs into row 1
    std::sort(arrayHaplo_count.begin(), arrayHaplo_count.end(),
              [](const std::vector<int> &a, const std::vector<int> &b)
              {
                  return a.size() > b.size();
              });

    num_rows_aoa = arrayHaplo_count.size();

    matrixAOA1 = new int *[num_rows_aoa];
    for (int i = 0; i < num_rows_aoa; i++)
    {

        matrixAOA1[i] = new int[num_pos];
        for (int k = 0; k < num_pos; k++)
        {
            matrixAOA1[i][k] = matrixHaploData[arrayHaplo_count[i][0]][k];
        }
    
    }

    for (int t = 0; t < num_rows_aoa; t++)
    {

        std::vector<int> chr_num_3D2sub;
        chr_num_3D2sub.insert(chr_num_3D2sub.begin(), arrayHaplo_count[t].begin(), arrayHaplo_count[t].end());
        std::vector<std::vector<std::string>> ids_conversion = this->convert3DChr2Sub(chr_num_3D2sub);

        aoa_chrs.push_back(chr_num_3D2sub);
        aoa_ids.push_back(ids_conversion[0]);
        aoa_2chrs.push_back(ids_conversion[1]);
    }

}

/// converts chr# from 3D matrix into subject ids
std::vector<std::vector<std::string>> HaploShare::convert3DChr2Sub(std::vector<int> &chr_num_3D)
{

    std::sort(chr_num_3D.begin(), chr_num_3D.end());

    std::vector<std::vector<std::string>> result;
    std::vector<std::string> temp1;
    std::vector<std::string> temp2;

    for (size_t k = 0; k < chr_num_3D.size(); k++)
    {

        for (size_t i = 0; i < ped1->cur_ped.size(); i++)
        {

            if (ped1->cur_ped[i].mat_3D_pos == chr_num_3D[k] ||
                ped1->cur_ped[i].pat_3D_pos == chr_num_3D[k])
            {
                temp1.push_back(ped1->cur_ped[i].subject_id);
                goto LB1;
            }
        }

    LB1:;
    }

    std::sort(temp1.begin(), temp1.end());

    temp2 = this->getSubsWithBothChrs(temp1);

    temp1.erase(unique(temp1.begin(), temp1.end()), temp1.end());

    result.push_back(temp1);
    result.push_back(temp2);

    return result;
}

/// returns duplicate ids
std::vector<std::string> HaploShare::getSubsWithBothChrs(std::vector<string> &duplicateIds)
{

    std::vector<std::string> result;

    std::map<std::string, int> countIds;

    for (auto &el : duplicateIds)
    {
        auto temp1 = countIds.insert(std::pair<std::string, int>(el, 1));
        if (temp1.second == false) // map key already exists
        {
            temp1.first->second++;
        }
    }

    for (auto &el : countIds)
    {
        if (el.second > 1)
        {
            result.push_back(el.first);
        }
    }

    std::sort(result.begin(), result.end());

    return result;
}

/// AOA2 matrix where 1s (REF alleles) replaced by row # starting from 1
void HaploShare::makeAOA2()
{

    matrixAOA2 = new int *[num_rows_aoa];
    for (int i = 0; i < num_rows_aoa; i++)
    {

        matrixAOA2[i] = new int[num_pos];
        for (int k = 0; k < num_pos; k++)
        {

            if (matrixAOA1[i][k] == 2)
            {
                matrixAOA2[i][k] = 0;
            }

            else

            {
                matrixAOA2[i][k] = i + 1;
            }
           
        }
     
    }

}

/// AOA3 matrix where 2s (ALT alleles) replaced by row # starting from 1
void HaploShare::makeAOA3()
{

    matrixAOA3 = new int *[num_rows_aoa];
    for (int i = 0; i < num_rows_aoa; i++)
    {

        matrixAOA3[i] = new int[num_pos];
        for (int k = 0; k < num_pos; k++)
        {

            if (matrixAOA1[i][k] == 1)
            {
                matrixAOA3[i][k] = 0;
            }

            else

            {
                matrixAOA3[i][k] = i + 1;
            }
        }
   
    }

}


/// determines whether to start computations using matrix AOA2 or AOA3 which depends on whether WGS data has lots of ALT alleles or not
int HaploShare::selectRoute() // route 1 = AOA2 and route 2 = AOA3, also could be route 0
{

    int route_num = 0;

    int count_one = 0;
    int count_two = 0;

    for (int i = 0; i < num_rows_aoa; i++)
    {
        int temp_one = 0;
        int temp_two = 0;

        for (int j = 0; j < num_pos; j++)
        {
            if (matrixAOA1[i][j] == 1)
            {
                temp_one++;
            }

            else if (matrixAOA1[i][j] == 2)

            {
                temp_two++;
            }
        }

        if (temp_one >= temp_two)
        {
            count_one++;
        }

        else

        {
            count_two++;
        }
    }

    if (count_one >= count_two)
    {
        route_num = 1;
    }

    else

    {
        route_num = 2;
    }

    return route_num;
}

/// searches for haplotype sharing
void HaploShare::search4HaploShare(int route_num)
{

    this->initAOA4Subst(route_num);

    std::vector<int> rowNums4SharedHaplo1stCheck = this->checkHpSharingMatrixTemp1();

    if (rowNums4SharedHaplo1stCheck[0] != -1)
    {
        std::vector<int> noComb2test(1, -1);
        this->checkCandidateGrpMatch(rowNums4SharedHaplo1stCheck, noComb2test);
    }


    std::vector<int> colnum4comb = this->getColNumComb(route_num);

    Combinations colcombs(&colnum4comb);
    std::vector<std::vector<int>> *colcombs2test = colcombs.getAllCombs();

    for (size_t i = 0; i < colcombs2test->size(); i++)
    {

        std::vector<int> temp_col4comb;

        for (size_t j = 0; j < (*colcombs2test)[i].size(); j++)
        {
            temp_col4comb.push_back((*colcombs2test)[i][j]);
        }

        this->testColSubstAOA(route_num, temp_col4comb);
        
    }
}

/// creates temporary matrix AOA4 (copy of AOA2 or AOA3, depending on route #) to test column substitutions
void HaploShare::initAOA4Subst(int route_num)
{

    matrixTemp1 = new int *[num_rows_aoa];
    for (int i = 0; i < num_rows_aoa; i++)
    {

        matrixTemp1[i] = new int[num_pos];
        for (int k = 0; k < num_pos; k++)
        {

            if (route_num == 1)
            {
                matrixTemp1[i][k] = matrixAOA2[i][k];
            }

            else if (route_num == 2)
            {
                matrixTemp1[i][k] = matrixAOA3[i][k];
            }

        }
   
    }

}

/// determines which rows have no zeros (identifies whether there is haplotype sharing between any subjects)
std::vector<int> HaploShare::checkHpSharingMatrixTemp1()
{
    std::vector<int> result;

    for (int i = 0; i < num_rows_aoa; i++)
    {
        int count1 = 0;

        for (int k = 0; k < num_pos; k++)
        {

            if (matrixTemp1[i][k] == 0)
            {
                goto LB1;
            }
            else
            {
                count1++;
            }
        
        }

        if (count1 == num_pos)
        {
            result.push_back(i);
        }

    LB1:;
    }

    if (result.empty())
    {
        result.push_back(-1);
    }

    return result; // return of -1 means no haplotype sharing but otherwise will return row numbers with same haplotype sequences
}

/// performs initial check for haplotype sharing (checks whether haplotype sharing group includes all subjects from core set)
void HaploShare::checkCandidateGrpMatch(std::vector<int> &rowNums4SharedHaplo, std::vector<int> &comb2test)
{
    std::vector<std::string> ids2check;

    for (size_t i = 0; i < rowNums4SharedHaplo.size(); i++)
    {
        for (size_t j = 0; j < aoa_ids[rowNums4SharedHaplo[i]].size(); j++)
        {  
            ids2check.push_back(aoa_ids[rowNums4SharedHaplo[i]][j]); // getting ids for specified row numbers from AOA1 matrix
        }
    }



    Comp_arrays_venn<string> comp2candidate(&ids2check, &candidate_grp);
    std::vector<string> commonIds = comp2candidate.getIntersection(); // finds intersection between ids from candidate grp (file) and ids being evaluated


    if (commonIds.size() == candidate_grp.size()) // matched
    {

        std::vector<std::vector<int>> geno_consist4share = this->checkGeno4ProposedAl(rowNums4SharedHaplo);

        if (geno_consist4share[0][0] == 0) // genotypes are consistent with sharing pattern
        {

            std::vector<int> suggestedAls;

            for (size_t i = 0; i < geno_consist4share[1].size(); i++)
            {
                suggestedAls.push_back(geno_consist4share[1][i]);
            }

            std::vector<int> subContribOnly1Chr;

            for (size_t i = 0; i < geno_consist4share[2].size(); i++)
            {
                subContribOnly1Chr.push_back(geno_consist4share[2][i]);
            }

            this->setCurSharingInfo(comb2test, rowNums4SharedHaplo, suggestedAls, subContribOnly1Chr);
        }

    }
}

/// checks if proposed alleles consistent with genotype data (core set)
std::vector<std::vector<int>> HaploShare::checkGeno4ProposedAl(std::vector<int> &rowNums4SharedHaplo)
{

    std::vector<std::vector<int>> result;
    std::vector<int> var_alleles;
    std::vector<int> subContribOnly1Chr; // stores subject index for those who contribute 2 chrs to haplosharing but based on genotype can contribute only one chr

    std::vector<int> negative_outcome(1, 0); // if 0 than everything is OK but if 1 than failed

    std::vector<int> varPos2varify = this->getVarPos2VarifyAOA1(rowNums4SharedHaplo);

    for (size_t k = 0; k < varPos2varify.size(); k++)
    {

        if (varPos2varify[k] != -1)
        {

            if (varPos2varify[k] == 1)
            {

                for (size_t h = 0; h < candidate_grp.size(); h++)
                {
                    // no single 2/2 genotype for ind from candidate group
                    // all 1/1 genotype for ind from candidate group with two chrs contributed to shared haplotype
                    std::vector<std::string> cur_sub_geno = dense_genotypes->getGeno((*cur_window_bps)[k], ped1->getSubjectIndex(candidate_grp[h]));
                    if (cur_sub_geno[0].compare("2") == 0 && cur_sub_geno[1].compare("2") == 0)
                    {
                        negative_outcome[0] = 1;
                        goto LB1;
                    }

                    int num_chrs_shared = this->getNumChrsShared(candidate_grp[h]); // if 1 = two chrs; -1 = nothing

                    if (num_chrs_shared == 1)
                    {
                        if (cur_sub_geno[0].compare("1") != 0 && cur_sub_geno[1].compare("1") != 0)
                        {
                            subContribOnly1Chr.push_back(ped1->getSubjectIndex(candidate_grp[h]));
                        }
                    }
                }

                var_alleles.push_back(1);
            }

            else if (varPos2varify[k] == 2)
            {

                for (size_t h = 0; h < candidate_grp.size(); h++)
                {

                    // no single 1/1 genotype for ind from candidate group
                    // all 2/2 genotype for ind from candidate group with two chrs contributed to shared haplotype
                    std::vector<std::string> cur_sub_geno = dense_genotypes->getGeno((*cur_window_bps)[k], ped1->getSubjectIndex(candidate_grp[h]));
                    if (cur_sub_geno[0].compare("1") == 0 && cur_sub_geno[1].compare("1") == 0)
                    {
                        negative_outcome[0] = 1;
                        goto LB1;
                    }

                    int num_chrs_shared = this->getNumChrsShared(candidate_grp[h]); // if 1 = two chrs; -1 = nothing

                    if (num_chrs_shared == 1)
                    {
                        if (cur_sub_geno[0].compare("2") != 0 && cur_sub_geno[1].compare("2") != 0)
                        {
                            subContribOnly1Chr.push_back(ped1->getSubjectIndex(candidate_grp[h]));
                        }
                    }
                }

                var_alleles.push_back(2);
            }

            else // case when variant has Zeros on haplotypes among all sharing Subjs
            {

                int proposed_al = 0;
                int count_ones_hm = 0;
                int count_twos_hm = 0;

                for (size_t h = 0; h < candidate_grp.size(); h++)
                {

                    int num_chrs_shared = this->getNumChrsShared(candidate_grp[h]); // if 1 = two chrs; -1 = nothing
                    std::vector<std::string> cur_sub_geno = dense_genotypes->getGeno((*cur_window_bps)[k], ped1->getSubjectIndex(candidate_grp[h]));

                    if (num_chrs_shared == 1)
                    {

                        if ((cur_sub_geno[0].compare("1") == 0 && cur_sub_geno[1].compare("2") == 0) || (cur_sub_geno[0].compare("2") == 0 && cur_sub_geno[1].compare("1") == 0))
                        {
                            subContribOnly1Chr.push_back(ped1->getSubjectIndex(candidate_grp[h]));
                        }
                    }

                    if (cur_sub_geno[0].compare("1") == 0 && cur_sub_geno[1].compare("1") == 0)
                    {
                        count_ones_hm++;
                    }

                    if (cur_sub_geno[0].compare("2") == 0 && cur_sub_geno[1].compare("2") == 0)
                    {
                        count_twos_hm++;
                    }
                }

                if (count_ones_hm > 0 && count_twos_hm > 0)
                {
                    negative_outcome[0] = 1;
                    goto LB1;
                }

                else if (count_ones_hm > 0)
                {
                    proposed_al = 1;
                }

                else if (count_twos_hm > 0)
                {
                    proposed_al = 2;
                }

                var_alleles.push_back(proposed_al);
            }
        }

        else

        {
            var_alleles.push_back(-1);
        }
    }

    LB1:; // way out of loop
    result.push_back(negative_outcome);

    result.push_back(var_alleles);

    std::sort(subContribOnly1Chr.begin(), subContribOnly1Chr.end());
    subContribOnly1Chr.erase(unique(subContribOnly1Chr.begin(), subContribOnly1Chr.end()), subContribOnly1Chr.end());
    result.push_back(subContribOnly1Chr); // was determined based on their genotype data


    return result; // V2[0], 1st el if = 0 then everything is OK, if = 1 then rejection; V2[1], 20 el are proposed alleles for 20 var positions if outcome is positive (== 0); V2[2], contains subject index for those who share only one haplotype with risk haplotype
}

/// determines at which variant positions on shared haplotype there are unresolved alleles in some subjects
std::vector<int> HaploShare::getVarPos2VarifyAOA1(std::vector<int> &rowNums4SharedHaplo)
{

    std::vector<int> result; // vector size equal to haplowindow size; -1 = nothing, 0 = all zeros in col, 1 = zeros and ones, 2 = zeros and twos

    for (int j = 0; j < num_pos; j++)
    {

        int count_zeros = 0;
        int count_ones = 0;
        int count_twos = 0;

        for (size_t i = 0; i < rowNums4SharedHaplo.size(); i++)
        {

            if (matrixAOA1[rowNums4SharedHaplo[i]][j] == 0)
            {
                count_zeros++;
            }

            else if (matrixAOA1[rowNums4SharedHaplo[i]][j] == 1)
            {
                count_ones++;
            }

            else if (matrixAOA1[rowNums4SharedHaplo[i]][j] == 2)
            {
                count_twos++;
            }
        }

        if (count_zeros > 0)
        {

            if (count_ones > 0)
            {
                result.push_back(1);
            }

            else if (count_twos > 0)
            {
                result.push_back(2);
            }

            else
            {
                result.push_back(0);
            }
        }

        else

        {
            result.push_back(-1);
        }
    }

    return result;
}

/// checks whether subject shares two chrs
int HaploShare::getNumChrsShared(std::string &sub_id)
{

    int result = -1;

    for (size_t i = 0; i < aoa_2chrs.size(); i++)
    {

        for (size_t j = 0; j < aoa_2chrs[i].size(); j++)
        {

            if (sub_id.compare(aoa_2chrs[i][j]) == 0)
            {
                result = 1;
                goto LB1;
            }
        }
    }

LB1:;

    return result;
}

/// recording current haplotype sharing results
void HaploShare::setCurSharingInfo(std::vector<int> &subsetCols,
std::vector<int> &sharingRows, std::vector<int> &suggestedAls,
std::vector<int> &subContribOnly1Chr)
{

    std::vector<std::vector<int>> temp1(6);

    temp1[0].insert(temp1[0].end(), subsetCols.begin(), subsetCols.end()); // combination of columns used to test substitution or -1 (none as during initial check)
    temp1[1].insert(temp1[1].end(), sharingRows.begin(), sharingRows.end()); // AOA1 matrix rows associated with shared haplotype
    temp1[2].insert(temp1[2].end(), suggestedAls.begin(), suggestedAls.end()); // suggested shared haplotype alleles
    temp1[3].insert(temp1[3].end(), subContribOnly1Chr.begin(), subContribOnly1Chr.end()); // determined based on genotype data


    cur_sharing.push_back(temp1);
}

/// identifies column # (variant position #) where substitutions will be performed to identify max # of chrs sharing same haplotype
std::vector<int> HaploShare::getColNumComb(int route_num)
{

    std::vector<int> result1;

    for (int k = 0; k < num_pos; k++)
    {

        int count_zeros = 0;
        int count_nonzeros = 0;

        for (int i = 0; i < num_rows_aoa; i++)
        {

            if (route_num == 1) // use AOA2
            {

                if (matrixAOA2[i][k] == 0)
                {
                    count_zeros++;
                }

                if (matrixAOA3[i][k] != 0)
                {
                    count_nonzeros++;
                }
            }

            else if (route_num == 2) // use AOA3
            {

                if (matrixAOA3[i][k] == 0)
                {
                    count_zeros++;
                }

                if (matrixAOA2[i][k] != 0)
                {
                    count_nonzeros++;
                }
            }
        }

        if (count_zeros > 0 || count_nonzeros > 0)
        {
            result1.push_back(k);
        }
    }

    return result1;
}

/// tests a combination of column number substitutions 
void HaploShare::testColSubstAOA(int route_num, std::vector<int> &comb2test)
{

    for (size_t j = 0; j < comb2test.size(); j++)
    {
        for (int i = 0; i < num_rows_aoa; i++)
        {
            if (route_num == 1)
            {
                matrixTemp1[i][comb2test[j]] = matrixAOA3[i][comb2test[j]];
            }

            else if (route_num == 2)
            {
                matrixTemp1[i][comb2test[j]] = matrixAOA2[i][comb2test[j]];
            }
        }
    }

    // Evaluate substitution
    std::vector<int> rowNums4SharedHaplo = this->checkHpSharingMatrixTemp1();

    if (rowNums4SharedHaplo[0] != -1)
    {

        this->checkCandidateGrpMatch(rowNums4SharedHaplo, comb2test);

    }
    // End of substitution evaluation

    // resetting Temp1 matrix used to test column substitutions
    for (size_t j = 0; j < comb2test.size(); j++)
    {
        for (int i = 0; i < num_rows_aoa; i++)
        {

            if (route_num == 1)
            {
                matrixTemp1[i][comb2test[j]] = matrixAOA2[i][comb2test[j]];
            }

            else if (route_num == 2)
            {
                matrixTemp1[i][comb2test[j]] = matrixAOA3[i][comb2test[j]];
            }
        }
    }

}

/// selects best haplotype sharing for a genomic window
void HaploShare::selectBestCurHaploShare() // execute only if cur_sharing is not empty!!!
{

    if (!cur_sharing.empty())
    {
        std::vector<int> els_4_step2 = this->selectCurHaploSharesStep1();

        std::vector<int> els_4_step3 = this->selectCurHaploSharesStep2(els_4_step2);

        std::vector<int> els_4_step4 = this->selectCurHaploSharesStep3(els_4_step3);

        this->selectCurHaploSharesStep4(els_4_step4);
    }
}

std::vector<int> HaploShare::selectCurHaploSharesStep1() // execute only if cur_sharing is not empty!!!
// purpose to get cases with the highest number of -1 and the lowest number of 0 at variant positions within genomic window
{

    std::vector<int> all_m_ones; // all -1
    std::vector<int> all_zeros;

    for (size_t i = 0; i < cur_sharing.size(); i++)
    {

        int count_m_ones = 0;
        int count_zeros = 0;

        for (size_t k = 0; k < cur_sharing[i][2].size(); k++)
        {

            if (cur_sharing[i][2][k] == -1)
            {
                count_m_ones++;
            }

            else if (cur_sharing[i][2][k] == 0)
            {
                count_zeros++;
            }
        }

        all_m_ones.push_back(count_m_ones);
        all_zeros.push_back(count_zeros);
    }

    int max_m_ones;
    std::vector<int> max_m_ones_els;
    std::vector<int> zeros_4_max_m_ones;

    if (!all_m_ones.empty())
    {
        max_m_ones = *std::max_element(all_m_ones.begin(), all_m_ones.end());

        for (size_t i = 0; i < all_m_ones.size(); i++)
        {
            if (all_m_ones[i] == max_m_ones)
            {
                max_m_ones_els.push_back(i);
                zeros_4_max_m_ones.push_back(all_zeros[i]);
            }
        }
    }

    int min_zeros;
    std::vector<int> els_2_next_step;

    if (!zeros_4_max_m_ones.empty())
    {
        min_zeros = *std::min_element(zeros_4_max_m_ones.begin(), zeros_4_max_m_ones.end());

        for (size_t i = 0; i < max_m_ones_els.size(); i++)
        {
            if (all_zeros[max_m_ones_els[i]] == min_zeros)
            {
                els_2_next_step.push_back(max_m_ones_els[i]);
            }
        }
    }

    return els_2_next_step;
}

std::vector<int> HaploShare::selectCurHaploSharesStep2(std::vector<int> &els_4_step2) // this function selects only unique haplotype sharings (unique row combinations, unique groups of haplotypes)
{

    std::vector<int> els_4_step3; // this is result that will be used in next step 3
    std::vector<std::vector<int>> rows2test;

    std::vector<int> temp1;
    for (size_t i = 0; i < cur_sharing[els_4_step2[0]][1].size(); i++)
    {
        temp1.push_back(cur_sharing[els_4_step2[0]][1][i]);
    }
    rows2test.push_back(temp1);
    els_4_step3.push_back(els_4_step2[0]);

    for (size_t i = 0; i < els_4_step2.size(); i++)
    {

        std::vector<int> temp2;

        for (size_t j = 0; j < cur_sharing[els_4_step2[i]][1].size(); j++)
        {
            temp2.push_back(cur_sharing[els_4_step2[i]][1][j]);
        }

        int status1 = 0;// 0 = combination new; 1 = similar combination exists
        int size = rows2test.size();
        for (int k = 0; k < size; k++)
        {
            std::vector<int> temp3;

            for (size_t h = 0; h < rows2test[k].size(); h++)
            {
                temp3.push_back(rows2test[k][h]);
            }

            Comp_arrays_venn<int> cur_comp(&temp3, &temp2);
            std::vector<int> comp_res = cur_comp.getDifference();

            if (comp_res.size() == 0)
            {
                status1 = 1;
            }
        }

        if (status1 == 0)
        {
            rows2test.push_back(temp2);
            els_4_step3.push_back(els_4_step2[i]);
        }

    }

    return els_4_step3;
}

/// selects entries where core group of cases shares haplotype with largest number of other subjects
std::vector<int> HaploShare::selectCurHaploSharesStep3(std::vector<int> &els_4_step3) // els for V1 of cur_sharing multi-D vector
{

    std::vector<int> numOthersSharing;
    std::vector<int> numOthersSharing_el;

    for (size_t i = 0; i < els_4_step3.size(); i++)
    {

        std::vector<int> rowNums4SharedHaplo;

        for (size_t k = 0; k < cur_sharing[els_4_step3[i]][1].size(); k++)
        {
            rowNums4SharedHaplo.push_back(cur_sharing[els_4_step3[i]][1][k]);
        }

        std::vector<std::string> all_ids_contrib2ShHp;

        for (size_t h = 0; h < rowNums4SharedHaplo.size(); h++)
        {

            for (size_t f = 0; f < aoa_ids[rowNums4SharedHaplo[h]].size(); f++)
            {
                all_ids_contrib2ShHp.push_back(aoa_ids[rowNums4SharedHaplo[h]][f]);
            }
        }

    
        std::sort(all_ids_contrib2ShHp.begin(), all_ids_contrib2ShHp.end());

        Comp_arrays_venn<string> mytest(&candidate_grp, &all_ids_contrib2ShHp);
        std::vector<string> extra_ids_who_share = mytest.getDifference();

        if (!extra_ids_who_share.empty())
        {

            std::vector<std::vector<int>> extraIdsShareConf = this->checkGeno4ProposedAlOthers(extra_ids_who_share, els_4_step3[i]);
            
            if (extraIdsShareConf[0][0] == -1)
            {
                numOthersSharing.push_back(0);
                numOthersSharing_el.push_back(i);
            }
            else
            {

                numOthersSharing.push_back(extraIdsShareConf[0].size());
                numOthersSharing_el.push_back(i);

                std::vector<int> sharingIdsOut;
                for (size_t h = 0; h < extraIdsShareConf[0].size(); h++)
                {
                    sharingIdsOut.push_back(extraIdsShareConf[0][h]);
                }
                cur_sharing[els_4_step3[i]][4].insert(cur_sharing[els_4_step3[i]][4].end(), sharingIdsOut.begin(), sharingIdsOut.end());

                std::vector<int> idsContrib1Chr;
                for (size_t h = 0; h < extraIdsShareConf[1].size(); h++)
                {
                    idsContrib1Chr.push_back(extraIdsShareConf[1][h]);
                }
                cur_sharing[els_4_step3[i]][5].insert(cur_sharing[els_4_step3[i]][5].end(), idsContrib1Chr.begin(), idsContrib1Chr.end());
            }
        }
    }

    int maxOthersShare;
    std::vector<int> el4maxOthersShare{-1};

    if (!numOthersSharing.empty())
    {
        maxOthersShare = *std::max_element(numOthersSharing.begin(), numOthersSharing.end());

        for (size_t i = 0; i < numOthersSharing.size(); i++)
        {
            if (numOthersSharing[i] == maxOthersShare)
            {
                el4maxOthersShare[0] = els_4_step3[numOthersSharing_el[i]];
                goto LB1;
            }
        }
    LB1:;
    }

    if (el4maxOthersShare[0] == -1)
    {
        el4maxOthersShare.push_back(els_4_step3[0]);
    }

    return el4maxOthersShare; // [0] = -1 if no other ids share haplotype with core group of cases or el # to use if there is a sharing; [1] = el # to use if [0] = -1
}

/// check if proposed alleles consistent with genotype data (core set and others)
std::vector<std::vector<int>> HaploShare::checkGeno4ProposedAlOthers(std::vector<std::string> &othersWhoShare, int curShareInst)
{

    std::vector<std::vector<int>> result;
    std::vector<int> subContribOnly1Chr; // stores subject index of those who contribute only one chr to shared haplotype

    std::vector<int> SharedHaploAls;
    for (size_t k = 0; k < cur_sharing[curShareInst][2].size(); k++)
    {
        SharedHaploAls.push_back(cur_sharing[curShareInst][2][k]);
    }

    for (size_t k = 0; k < SharedHaploAls.size(); k++)
    {

        if (SharedHaploAls[k] != -1)
        {

            if (SharedHaploAls[k] == 1)
            {
                std::vector<int> el2remove;

                for (size_t h = 0; h < othersWhoShare.size(); h++)
                {
                    // no single 2/2 genotype for ind from others group
                    // all 1/1 genotype for ind from others group with two chrs contributed to shared haplotype
                    int num_chrs_shared;
                    std::vector<std::string> cur_sub_geno = dense_genotypes->getGeno((*cur_window_bps)[k], ped1->getSubjectIndex(othersWhoShare[h]));
                    if (cur_sub_geno[0].compare("2") == 0 && cur_sub_geno[1].compare("2") == 0)
                    {
                        el2remove.push_back(h);
                        goto LB1;
                    }

                    num_chrs_shared = this->getNumChrsShared(othersWhoShare[h]); // if 1 = two chrs; -1 = nothing

                    if (num_chrs_shared == 1)
                    {
                        if (cur_sub_geno[0].compare("1") != 0 && cur_sub_geno[1].compare("1") != 0)
                        {
                            subContribOnly1Chr.push_back(ped1->getSubjectIndex(othersWhoShare[h]));
                        }
                    }

                LB1:;
                }

                if (!el2remove.empty())
                {

                    std::sort(el2remove.begin(), el2remove.end());
                    for (int t = el2remove.size() - 1; t >= 0; t--)
                    {
                        othersWhoShare.erase(othersWhoShare.begin() + el2remove[t]);
                    }
                }

            }

            else if (SharedHaploAls[k] == 2)
            {
                
                std::vector<int> el2remove;

                for (size_t h = 0; h < othersWhoShare.size(); h++)
                {

                    // no single 1/1 genotype for ind from others group
                    // all 2/2 genotype for ind from others group with two chrs contributed to shared haplotype
                    int num_chrs_shared;
                    std::vector<std::string> cur_sub_geno = dense_genotypes->getGeno((*cur_window_bps)[k], ped1->getSubjectIndex(othersWhoShare[h]));
                    if (cur_sub_geno[0].compare("1") == 0 && cur_sub_geno[1].compare("1") == 0)
                    {
                        el2remove.push_back(h);
                        goto LB2;
                    }

                    num_chrs_shared = this->getNumChrsShared(othersWhoShare[h]); // if 1 = two chrs; -1 = nothing

                    if (num_chrs_shared == 1)
                    {
                        if (cur_sub_geno[0].compare("2") != 0 && cur_sub_geno[1].compare("2") != 0)
                        {
                            subContribOnly1Chr.push_back(ped1->getSubjectIndex(othersWhoShare[h]));
                        }
                    }

                LB2:;
                }


                if (!el2remove.empty())
                {

                    std::sort(el2remove.begin(), el2remove.end());
                    for (int t = el2remove.size() - 1; t >= 0; t--)
                    {
                        othersWhoShare.erase(othersWhoShare.begin() + el2remove[t]);
                    }
                }

            }

            else if (SharedHaploAls[k] == 0) // case when variant has Zeros on haplotypes among all sharing Subjs
            {

                for (size_t h = 0; h < othersWhoShare.size(); h++)
                {

                    int num_chrs_shared = this->getNumChrsShared(othersWhoShare[h]); // if 1 = two chrs; -1 = nothing
                    std::vector<std::string> cur_sub_geno = dense_genotypes->getGeno((*cur_window_bps)[k], ped1->getSubjectIndex(othersWhoShare[h]));

                    if (num_chrs_shared == 1)
                    {

                        if ((cur_sub_geno[0].compare("1") == 0 && cur_sub_geno[1].compare("2") == 0) || (cur_sub_geno[0].compare("2") == 0 && cur_sub_geno[1].compare("1") == 0))
                        {
                            subContribOnly1Chr.push_back(ped1->getSubjectIndex(othersWhoShare[h]));
                        }
                    }
                }
            }
        }
    }


    if (!othersWhoShare.empty())
    {

        std::vector<int> temp1;
        for (size_t g = 0; g < othersWhoShare.size(); g++)
        {
            int temp2 = ped1->getSubjectIndex(othersWhoShare[g]);
            temp1.push_back(temp2);
        }

        result.push_back(temp1);
    }
    else
    {
        std::vector<int> temp1(1, -1);
        result.push_back(temp1);
    }

    if (result[0][0] != -1)
    {
        std::sort(subContribOnly1Chr.begin(), subContribOnly1Chr.end());
        subContribOnly1Chr.erase(unique(subContribOnly1Chr.begin(), subContribOnly1Chr.end()), subContribOnly1Chr.end());
        result.push_back(subContribOnly1Chr);
    }

    return result; // V2[0], 1st el if = -1 then rejection, otherwise = ids of others who share; V2[1], contains subject index for those not from core group of cases who share only one haplotype with risk haplotype
}

/// mostly recording resutls of haplotype sharing for a genomic window
void HaploShare::selectCurHaploSharesStep4(std::vector<int> &els_4_step4)
{

    int el_4_step4;

    if (els_4_step4[0] != -1)// we have others who share haplotype with cases from core group
    {
        el_4_step4 = els_4_step4[0];
    }
    else
    {
        el_4_step4 = els_4_step4[1];// pushing 1st element from step2
    }

    std::vector<std::vector<int>> result1(4);

    std::vector<int> sharedHaploRows;
    for (size_t i = 0; i < cur_sharing[el_4_step4][1].size(); i++)
    {
        sharedHaploRows.push_back(cur_sharing[el_4_step4][1][i]);
    }

    std::vector<int> shared_chr_3D;
    for (size_t i = 0; i < sharedHaploRows.size(); i++)
    {

        for (size_t g = 0; g < aoa_chrs[sharedHaploRows[i]].size(); g++)
        {
            shared_chr_3D.push_back(aoa_chrs[sharedHaploRows[i]][g]);
        }
    }

    std::vector<int> shared_chr_3D_CandGrp;
    for (size_t i = 0; i < candidate_grp.size(); i++)
    {

        for (size_t j = 0; j < shared_chr_3D.size(); j++)
        {

            if (shared_chr_3D[j] == ped1->cur_ped[ped1->getSubjectIndex(candidate_grp[i])].mat_3D_pos)
            {
                shared_chr_3D_CandGrp.push_back(shared_chr_3D[j]);
            }

            else if (shared_chr_3D[j] == ped1->cur_ped[ped1->getSubjectIndex(candidate_grp[i])].pat_3D_pos)
            {
                shared_chr_3D_CandGrp.push_back(shared_chr_3D[j]);
            }
           
        }
     
    }

    std::vector<int> others;
    std::vector<int> shared_chr_3D_Others;

    if (els_4_step4[0] != -1)
    {
        
        if (!cur_sharing[el_4_step4][4].empty())
        {

            for (size_t i = 0; i < cur_sharing[el_4_step4][4].size(); i++)
            {

                others.push_back(cur_sharing[el_4_step4][4][i]);
            }
        }

        for (size_t i = 0; i < others.size(); i++)
        {

            for (size_t j = 0; j < shared_chr_3D.size(); j++)
            {

                if (shared_chr_3D[j] == ped1->cur_ped[others[i]].mat_3D_pos)
                {
                    shared_chr_3D_Others.push_back(shared_chr_3D[j]);
                }

                else if (shared_chr_3D[j] == ped1->cur_ped[others[i]].pat_3D_pos)
                {
                    shared_chr_3D_Others.push_back(shared_chr_3D[j]);
                }
            }
        }
    }

    std::vector<int> chr3D4contribOnly1ChrCandGrp;

    if (!cur_sharing[el_4_step4][3].empty())
    {

        for (size_t i = 0; i < cur_sharing[el_4_step4][3].size(); i++)
        {
            chr3D4contribOnly1ChrCandGrp.push_back(ped1->cur_ped[cur_sharing[el_4_step4][3][i]].mat_3D_pos);
            chr3D4contribOnly1ChrCandGrp.push_back(ped1->cur_ped[cur_sharing[el_4_step4][3][i]].pat_3D_pos);
        }
    }

    std::vector<int> chr3D4contribOnly1ChrOthers;

    if (els_4_step4[0] != -1)
    {

        if (!cur_sharing[el_4_step4][5].empty())
        {

            for (size_t i = 0; i < cur_sharing[el_4_step4][5].size(); i++)
            {
                chr3D4contribOnly1ChrOthers.push_back(ped1->cur_ped[cur_sharing[el_4_step4][5][i]].mat_3D_pos);
                chr3D4contribOnly1ChrOthers.push_back(ped1->cur_ped[cur_sharing[el_4_step4][5][i]].pat_3D_pos);
            }
        }
    }


    result1[0].insert(result1[0].end(), shared_chr_3D_CandGrp.begin(), shared_chr_3D_CandGrp.end());
    result1[1].insert(result1[1].end(), chr3D4contribOnly1ChrCandGrp.begin(), chr3D4contribOnly1ChrCandGrp.end());
    result1[2].insert(result1[2].end(), shared_chr_3D_Others.begin(), shared_chr_3D_Others.end());
    result1[3].insert(result1[3].end(), chr3D4contribOnly1ChrOthers.begin(), chr3D4contribOnly1ChrOthers.end());

    global_hp_share->setValue(cur_wind, result1);

    std::vector<int> sharedHaploTemplate;
    for (size_t i = 0; i < cur_sharing[el_4_step4][2].size(); i++)
    {
        sharedHaploTemplate.push_back(cur_sharing[el_4_step4][2][i]);
    }

    for (size_t j = 0; j < sharedHaploTemplate.size(); j++)
    {
        if (sharedHaploTemplate[j] == -1)
        {
            sharedHaploOut.push_back(matrixAOA1[sharedHaploRows[0]][j]);
        }

        else

        {
            sharedHaploOut.push_back(sharedHaploTemplate[j]);
        }
    }

    *sharedHaploStream << chr << ":" << (*cur_window_bps)[0] << "-" << (*cur_window_bps)[cur_window_bps->size() - 1] << " ";
    for (size_t j = 0; j < sharedHaploOut.size(); j++)
    {
        *sharedHaploStream << sharedHaploOut[j];
    }
    *sharedHaploStream << std::endl;
   
}

/// possibility for a second shared haplotype determined based on genotype data in subjects from core set
int HaploShare::check2ndSharedHaploPossibility()
{
    int result = -1;

    int het_pos_num = 0;
    int hm1_pos_num = 0;
    int hm2_pos_num = 0;
    int miss_pos_num = 0;
    int genos_total = candidate_grp.size();
    int geno_pos_sum = 0;

    for (int k = 0; k < num_pos; k++)
    {

        int count_g11 = 0;
        int count_g22 = 0;
        int count_g12 = 0;
        int count_miss1and1 = 0;
        int count_miss1and2 = 0;
        int count_miss2 = 0;

        for (size_t h = 0; h < candidate_grp.size(); h++)
        {

            std::vector<std::string> cur_sub_geno = dense_genotypes->getGeno((*cur_window_bps)[k], ped1->getSubjectIndex(candidate_grp[h]));

            if (cur_sub_geno[0].compare("1") == 0 && cur_sub_geno[1].compare("1") == 0)
            {
                count_g11++;
            }

            else if (cur_sub_geno[0].compare("2") == 0 && cur_sub_geno[1].compare("2") == 0)
            {
                count_g22++;
            }

            else if ((cur_sub_geno[0].compare("1") == 0 && cur_sub_geno[1].compare("2") == 0) ||
                     (cur_sub_geno[0].compare("2") == 0 && cur_sub_geno[1].compare("1") == 0))
            {
                count_g12++;
            }

            else if (cur_sub_geno[0].compare("0") == 0 && cur_sub_geno[1].compare("0") == 0)
            {
                count_miss2++;
            }

            else if ((cur_sub_geno[0].compare("0") == 0 && cur_sub_geno[1].compare("1") == 0) ||
                     (cur_sub_geno[1].compare("0") == 0 && cur_sub_geno[0].compare("1") == 0))
            {
                count_miss1and1++;
            }

            else if ((cur_sub_geno[0].compare("0") == 0 && cur_sub_geno[1].compare("2") == 0) ||
                     (cur_sub_geno[1].compare("0") == 0 && cur_sub_geno[0].compare("2") == 0))
            {
                count_miss1and2++;
            }
        }

        if (count_g11 == 0 && count_g12 == 0 && count_miss1and1 == 0 && count_miss2 < genos_total)
        {
            hm2_pos_num++;
        }

        else if (count_g22 == 0 && count_g12 == 0 && count_miss1and2 == 0 && count_miss2 < genos_total)
        {
            hm1_pos_num++;
        }

        else if (count_g11 == 0 && count_g22 == 0 && count_miss2 < genos_total)
        {
            het_pos_num++;
        }

        else if (count_miss2 == genos_total)
        {
            miss_pos_num++;
        }
    }

    geno_pos_sum = het_pos_num + hm1_pos_num + hm2_pos_num + miss_pos_num;

    if (het_pos_num > 0 && geno_pos_sum == num_pos)
    {
        result = 1;
    }

    return result;
}

/// determines sequence of the second shared haplotype using sequence of the first shared haplotype and genotype data
void HaploShare::determine2ndSharedHaploSeq()
{

    int genos_total = candidate_grp.size(); // # of genotypes at a variant position in core set (equals # of sbujects in core set)
    std::vector<std::string> second_shared_haplo;
    std::vector<int> second_shared_haplo_int;
    int check1 = 0;// counts how many variant positions will have different alleles on two shared haplotypes
    
    for (int k = 0; k < num_pos; k++)
    {

        int cur_al_hp1 = sharedHaploOut[k];
        int cur_al_hp2 = 0;

        int count_g11 = 0;
        int count_g22 = 0;
        int count_g12 = 0;
        int count_miss1and1 = 0;
        int count_miss1and2 = 0;
        int count_miss2 = 0;

        for (size_t h = 0; h < candidate_grp.size(); h++)
        {

            std::vector<std::string> cur_sub_geno = dense_genotypes->getGeno((*cur_window_bps)[k], ped1->getSubjectIndex(candidate_grp[h]));

            if (cur_sub_geno[0].compare("1") == 0 && cur_sub_geno[1].compare("1") == 0)
            {
                count_g11++;
            }

            else if (cur_sub_geno[0].compare("2") == 0 && cur_sub_geno[1].compare("2") == 0)
            {
                count_g22++;
            }

            else if ((cur_sub_geno[0].compare("1") == 0 && cur_sub_geno[1].compare("2") == 0) ||
                     (cur_sub_geno[0].compare("2") == 0 && cur_sub_geno[1].compare("1") == 0))
            {
                count_g12++;
            }

            else if (cur_sub_geno[0].compare("0") == 0 && cur_sub_geno[1].compare("0") == 0)
            {
                count_miss2++;
            }

            else if ((cur_sub_geno[0].compare("0") == 0 && cur_sub_geno[1].compare("1") == 0) ||
                     (cur_sub_geno[1].compare("0") == 0 && cur_sub_geno[0].compare("1") == 0))
            {
                count_miss1and1++;
            }

            else if ((cur_sub_geno[0].compare("0") == 0 && cur_sub_geno[1].compare("2") == 0) ||
                     (cur_sub_geno[1].compare("0") == 0 && cur_sub_geno[0].compare("2") == 0))
            {
                count_miss1and2++;
            }
        }

        if (count_g11 == 0 && count_g12 == 0 && count_miss1and1 == 0 && count_miss2 < genos_total && count_g22 > 0)
        {
            if (cur_al_hp1 == 2)
            {
                cur_al_hp2 = 2;
            }
        }

        else if (count_g22 == 0 && count_g12 == 0 && count_miss1and2 == 0 && count_miss2 < genos_total && count_g11 > 0)
        {
            if (cur_al_hp1 == 1)
            {
                cur_al_hp2 = 1;
            }
        }
        
        else if (count_g11 == 0 && count_g22 == 0 && count_miss2 < genos_total && (count_g12 > 0 || count_miss1and1 > 0 || count_miss1and2 > 0))
        {
            check1++;
            
            if (cur_al_hp1 == 1)
            {
                cur_al_hp2 = 2;
            }

            else if (cur_al_hp1 == 2)
            {
                cur_al_hp2 = 1;
            }
        }

        else if (count_miss2 == genos_total)
        {
            cur_al_hp2 = 0;
        }
        
        else if (count_miss2 + count_miss1and1 == genos_total && count_miss1and2 == 0)
        {
            if (cur_al_hp1 == 1)
            {
                cur_al_hp2 = 0;
            }
            else
            {
                cur_al_hp2 = 1;
            }
            
        }
        
        else if (count_miss2 + count_miss1and2 == genos_total && count_miss1and1 == 0)
        {
            if (cur_al_hp1 == 2)
            {
                cur_al_hp2 = 0;
            }
            else
            {
                cur_al_hp2 = 2;
            }
           
        }

        std::string cur_al_hp2_out = std::to_string(cur_al_hp2);
        second_shared_haplo.push_back(cur_al_hp2_out);
        second_shared_haplo_int.push_back(cur_al_hp2);
    }


    if (check1 > 0)// ensures that at least one var pos has het genotypes in all cases from core group
    {

        std::string cur_wind_el_str = std::to_string(cur_wind - 1); // window# - 1

        std::string windBpBounds = std::to_string(chr) + ":" + std::to_string((*cur_window_bps)[0]) + "-" + std::to_string((*cur_window_bps)[cur_window_bps->size() - 1]);

        global_hp_share->record2ndSharedHaplo(cur_wind_el_str, windBpBounds, second_shared_haplo);

        std::vector<int> candidate_grp3Dchrs = this->get3Dchrs4candidate_grp();

        global_hp_share->setValue2(cur_wind, candidate_grp3Dchrs);

        std::vector<int> othersSharing2ndHaplo = this->determineOthersWhoShare2ndHaplo(second_shared_haplo_int);
        
        global_hp_share->setValue3(cur_wind, othersSharing2ndHaplo);
    
    }
   
}


std::vector<int> HaploShare::get3Dchrs4candidate_grp()
{

    std::vector<int> result;

    for (size_t i = 0; i < candidate_grp.size(); i++)
    {

        result.push_back (ped1->cur_ped[ped1->getSubjectIndex(candidate_grp[i])].mat_3D_pos);

        result.push_back (ped1->cur_ped[ped1->getSubjectIndex(candidate_grp[i])].pat_3D_pos);

    }

    return result;

}


/// compares sequences of the second shared haplotype with subject haplotypes to identify others (non core set members) who share the haplotype
std::vector<int> HaploShare::determineOthersWhoShare2ndHaplo(std::vector<int> &second_shared_haplo_int)
{

    std::vector<int> othersSharing2ndHaplo; // contains 3D chrs

    std::vector<std::string> wgsDataSubjects;

    for (size_t i = 0; i < ped1->cur_ped.size(); ++i)
    {
        if (ped1->cur_ped[i].wgs_data == 1)
        {
            wgsDataSubjects.push_back(ped1->cur_ped[i].subject_id);
        }
    }

    Comp_arrays_venn<string> mytest(&candidate_grp, &wgsDataSubjects);
    std::vector<std::string> wgsDataSubOthers = mytest.getDifference();

    for (size_t i = 0; i < wgsDataSubOthers.size(); ++i)
    {

        int curSubInd = ped1->getSubjectIndex(wgsDataSubOthers[i]);
       
        std::vector<int> hp_share_status = global_hp_share->readHaploShareData(cur_wind, wgsDataSubOthers[i]); // [0] = maternal; [1] = paternal


        int success_count = 0;

        if (hp_share_status[0] == 3)
        {
            goto LB1;
        }

        
        for (int k = 0; k < num_pos; k++)
        {
            
            std::vector<std::string> cur_sub_geno = dense_genotypes->getGeno((*cur_window_bps)[k], curSubInd);

            if (second_shared_haplo_int[k] == 1 && (matrixHaploData[ped1->cur_ped[curSubInd].mat_3D_pos][k] == 1 ||
            matrixHaploData[ped1->cur_ped[curSubInd].mat_3D_pos][k] == 0) && ( 
            (cur_sub_geno[0].compare("1") == 0 && cur_sub_geno[1].compare("1") == 0) || 
            (cur_sub_geno[0].compare("1") == 0 && cur_sub_geno[1].compare("2") == 0) ||
            (cur_sub_geno[0].compare("2") == 0 && cur_sub_geno[1].compare("1") == 0) ||
            (cur_sub_geno[0].compare("0") == 0 && cur_sub_geno[1].compare("0") == 0)))
            {
                success_count++;
            }

            else if (second_shared_haplo_int[k] == 2 && (matrixHaploData[ped1->cur_ped[curSubInd].mat_3D_pos][k] == 2 ||
            matrixHaploData[ped1->cur_ped[curSubInd].mat_3D_pos][k] == 0) && (
            (cur_sub_geno[0].compare("2") == 0 && cur_sub_geno[1].compare("2") == 0) || 
            (cur_sub_geno[0].compare("1") == 0 && cur_sub_geno[1].compare("2") == 0) ||
            (cur_sub_geno[0].compare("2") == 0 && cur_sub_geno[1].compare("1") == 0) ||
            (cur_sub_geno[0].compare("0") == 0 && cur_sub_geno[1].compare("0") == 0)))
            {
                success_count++;
            }

            else if (second_shared_haplo_int[k] == 0 && (matrixHaploData[ped1->cur_ped[curSubInd].mat_3D_pos][k] == 0 || 
            matrixHaploData[ped1->cur_ped[curSubInd].mat_3D_pos][k] == 1 || matrixHaploData[ped1->cur_ped[curSubInd].mat_3D_pos][k] == 2))
            {
                success_count++;
            }
        }

        if (success_count == num_pos)
        {
            othersSharing2ndHaplo.push_back(ped1->cur_ped[curSubInd].mat_3D_pos);
        }

        LB1:;


        int success_count2 = 0;

        if (hp_share_status[1] == 3)
        {
            goto LB2;
        }

        

        for (int k = 0; k < num_pos; k++)
        {

            std::vector<std::string> cur_sub_geno = dense_genotypes->getGeno((*cur_window_bps)[k], curSubInd);

            if (second_shared_haplo_int[k] == 1 && (matrixHaploData[ped1->cur_ped[curSubInd].pat_3D_pos][k] == 1 ||
            matrixHaploData[ped1->cur_ped[curSubInd].pat_3D_pos][k] == 0) && (
            (cur_sub_geno[0].compare("1") == 0 && cur_sub_geno[1].compare("1") == 0) || 
            (cur_sub_geno[0].compare("1") == 0 && cur_sub_geno[1].compare("2") == 0) ||
            (cur_sub_geno[0].compare("2") == 0 && cur_sub_geno[1].compare("1") == 0) ||
            (cur_sub_geno[0].compare("0") == 0 && cur_sub_geno[1].compare("0") == 0)))
            {
                success_count2++;
            }

            else if (second_shared_haplo_int[k] == 2 && (matrixHaploData[ped1->cur_ped[curSubInd].pat_3D_pos][k] == 2 ||
            matrixHaploData[ped1->cur_ped[curSubInd].pat_3D_pos][k] == 0) && (
            (cur_sub_geno[0].compare("2") == 0 && cur_sub_geno[1].compare("2") == 0) || 
            (cur_sub_geno[0].compare("1") == 0 && cur_sub_geno[1].compare("2") == 0) ||
            (cur_sub_geno[0].compare("2") == 0 && cur_sub_geno[1].compare("1") == 0) ||
            (cur_sub_geno[0].compare("0") == 0 && cur_sub_geno[1].compare("0") == 0)))   
            {
                success_count2++;
            }

            else if (second_shared_haplo_int[k] == 0 && (matrixHaploData[ped1->cur_ped[curSubInd].pat_3D_pos][k] == 0 ||
            matrixHaploData[ped1->cur_ped[curSubInd].pat_3D_pos][k] == 1 || matrixHaploData[ped1->cur_ped[curSubInd].pat_3D_pos][k] == 2))
            {
                success_count2++;
            }
        }

        if (success_count2 == num_pos)
        {
            othersSharing2ndHaplo.push_back(ped1->cur_ped[curSubInd].pat_3D_pos);
        }
    
    LB2:;
    }

    return othersSharing2ndHaplo;
}

/// checks if there is a possibility for 2nd ambiguous shared haplotype 
std::vector<int> HaploShare::detectAmbigSharedHaplosPresence()
{
   
    std::vector<int> result1 = {-1};// -1 = nothing; number > 0 = possible ambiguouse haplotype sharing

    std::vector<int> chrs3D4SecondSharedHp;
    int count1 = 0;

    for (size_t i = 0; i < candidate_grp.size(); i++)
    {
        
        std::vector<int> share4CurSub = global_hp_share->readHaploShareData(cur_wind, candidate_grp[i]);
    
        if (share4CurSub[0] == -1 || share4CurSub[1] == -1 ||
        (share4CurSub[0] == 2 && share4CurSub[1] == 2))
        {
            count1++;

            if (share4CurSub[0] == -1)
            {
                chrs3D4SecondSharedHp.push_back(ped1->cur_ped[ped1->getSubjectIndex(candidate_grp[i])].mat_3D_pos);
            }

            else if (share4CurSub[1] == -1)
            {
                chrs3D4SecondSharedHp.push_back(ped1->cur_ped[ped1->getSubjectIndex(candidate_grp[i])].pat_3D_pos);
            }

            else if (share4CurSub[0] == 2 && share4CurSub[1] == 2)
            {
                chrs3D4SecondSharedHp.push_back(ped1->cur_ped[ped1->getSubjectIndex(candidate_grp[i])].mat_3D_pos);
                chrs3D4SecondSharedHp.push_back(ped1->cur_ped[ped1->getSubjectIndex(candidate_grp[i])].pat_3D_pos);
            }

        }
    
    }

    if (count1 == static_cast<int>(candidate_grp.size()))
    {
        return chrs3D4SecondSharedHp;

    }
    else
    {
        return result1;
    }

}


void HaploShare::evalAmbigSharedHaplosPresence(std::vector<int> &chrs3D4SecondSharedHp)
{

    std::vector<string> proposed_al4ambigHaplo;
    std::vector<int> proposed_al4ambigHaplo_int;
   
    for (int j = 0; j < num_pos; j++)
    {

        int count_0 = 0;
        int count_1 = 0;
        int count_2 = 0;
        
        for (size_t i = 0; i < chrs3D4SecondSharedHp.size(); i++)
        {

            if (matrixHaploData[chrs3D4SecondSharedHp[i]][j] == 0)
            {
                count_0++;
            }

            else if (matrixHaploData[chrs3D4SecondSharedHp[i]][j] == 1)
            {
                count_1++;
            }

            else if (matrixHaploData[chrs3D4SecondSharedHp[i]][j] == 2)
            {
                count_2++;
            }

        }


        if (count_1 > 0 && count_2 > 0)
        {
            goto LB1;
        }
        else
        {

            if (count_1 > 0)
            {
                proposed_al4ambigHaplo.push_back("1");
                proposed_al4ambigHaplo_int.push_back(1);
            }

            else if (count_2 > 0)
            {
                proposed_al4ambigHaplo.push_back("2");
                proposed_al4ambigHaplo_int.push_back(2);
            }

            else if (count_0 == static_cast<int>(chrs3D4SecondSharedHp.size()))
            {
                proposed_al4ambigHaplo.push_back("0");
                proposed_al4ambigHaplo_int.push_back(0);
            }

        }

    }


    if (proposed_al4ambigHaplo.size() == static_cast<size_t>(num_pos))
    {

        std::string cur_wind_el_str = std::to_string(cur_wind - 1); // window# - 1

        std::string windBpBounds = std::to_string(chr) + ":" + std::to_string((*cur_window_bps)[0]) + "-" + std::to_string((*cur_window_bps)[cur_window_bps->size() - 1]);

        global_hp_share->record2ndSharedHaplo(cur_wind_el_str, windBpBounds, proposed_al4ambigHaplo);

        global_hp_share->setValue2(cur_wind, chrs3D4SecondSharedHp);

        std::vector<int> othersSharing2ndHaplo = this->determineOthersWhoShare2ndHaplo(proposed_al4ambigHaplo_int);
        
        global_hp_share->setValue3(cur_wind, othersSharing2ndHaplo);

    }

    LB1:;

}


/// input functions
std::vector<int> split(const std::string &s, char delim)
{
    std::vector<int> elems;
    std::stringstream ss(s);
    std::string number;
    while (std::getline(ss, number, delim))
    {
        elems.push_back(std::stoi(number));
    }
    return elems;
}

std::vector<double> split3(const std::string &s, char delim)
{
    std::vector<double> elems;
    std::stringstream ss(s);
    std::string number;
    while (std::getline(ss, number, delim))
    {
        elems.push_back(std::stod(number));
    }
    return elems;
}

std::vector<string> split4(const std::string &s, char delim)
{
   
    std::vector<string> elems;
    std::stringstream ss(s);
    std::string number;
    while (std::getline(ss, number, delim))
    {
        number = std::regex_replace(number, std::regex("^ +| +$|( ) +"), "$1");
        elems.push_back(number);
    }
    return elems;
}

void trim(string& str)
{
   while(str[0] == ' ') str.erase(str.begin());
   while(str[str.size() - 1] == ' ') str.pop_back();
}

// Helper to trim whitespace from a string (both ends)
std::string trim(const std::string& s) {
    auto start = s.begin();
    while (start != s.end() && std::isspace(*start)) start++;

    auto end = s.end();
    do {
        end--;
    } while (std::distance(start, end) > 0 && std::isspace(*end));

    return std::string(start, end + 1);
}

void readParameterFile(std::string &parameterFileName, Parameters &current_pars)
{

    std::ifstream inParameterFile;
    inParameterFile.open(parameterFileName, std::ios::in);

    if (!inParameterFile)
    {
        cerr << "Parameter file could not be opened" << std::endl;
        exit(1);
    }

    
    char del1 = '#';
    std::string line;
 
    if (inParameterFile.is_open())
    {

        while (getline(inParameterFile, line))
        {
           
            trim(line);
            
            if (line.length() == 0)
            {
                continue;  
            }

            else if (line.at(0) == '*')
            {
                continue;  
            }

            else if(line.find("@") != line.npos)
            {
                continue;
            }

            if (line.find("#") != line.npos)
            {
                
                std::vector<string> result1 =  split4(line, del1);

                if ((result1[0].compare("1") == 0) && (result1.size() > 1)) // equal to each other
                {
                    current_pars.pedigree_file_path = result1[1];

                    if (result1[1].empty())
                    {
                        cerr << "Pedigree file path is not specified" << std::endl;
                        exit(1);
                    }
                }

                else if ((result1[0].compare("2") == 0) && (result1.size() > 1)) // equal to each other
                {
                    current_pars.dense_marker_pos_file_path = result1[1];
                }

                else if ((result1[0].compare("3") == 0) && (result1.size() > 1)) // equal to each other
                {
                    current_pars.dense_marker_geno_txt_file_path = result1[1];
                }

                else if ((result1[0].compare("4") == 0) && (result1.size() > 1)) // equal to each other
                {
                    trim(result1[1]);
                    regex rex1("\\s{2,}");
                    string del1 = " ";
                    string result2 = regex_replace(result1[1], rex1, del1);
                    std::vector<double> region = split3(result2, ' ');
                    current_pars.linkage_region = region; 
                }

                else if ((result1[0].compare("5") == 0) && (result1.size() > 1)) // equal to each other
                {
                    double result2 = std::stod(result1[1]);
                    current_pars.maxLODmarker = result2;
                }

                else if ((result1[0].compare("6") == 0) && (result1.size() > 1)) // equal to each other
                {
                    current_pars.output_dir_path = result1[1];
                }

                else if ((result1[0].compare("7") == 0) && (result1.size() > 1)) // equal to each other
                {
                    current_pars.framework_marker_pos_file_path = result1[1];
                }

                else if ((result1[0].compare("8") == 0) && (result1.size() > 1)) // equal to each other
                {
                    current_pars.meiosis_file_path = result1[1];
                }

                else if ((result1[0].compare("9") == 0) && (result1.size() > 1)) // equal to each other
                {
                    long result2 = std::stol(result1[1]);
                    current_pars.num_iter = result2;
                }

                else if ((result1[0].compare("10") == 0) && (result1.size() > 1)) // equal to each other
                {
                    current_pars.haplotype_file_path = result1[1];
                }

                else if ((result1[0].compare("11") == 0) && (result1.size() > 1)) // equal to each other
                {
                    current_pars.candidate_grp_file_path = result1[1];
                }

                else if ((result1[0].compare("12") == 0) && (result1.size() > 1)) // equal to each other
                {
                    long result2 = std::stol(result1[1]);
                    current_pars.seed = result2;
                }

            }
        
        }
       
       
        if (current_pars.run_type == 3)
        {
            current_pars.haplotype_file_path = current_pars.output_dir_path + "/haplotype_sequences.txt";
            current_pars.candidate_grp_file_path = current_pars.output_dir_path + "/core_set_of_cases.txt";
        }


    }
   
    inParameterFile.close();
}

void readPedigreeFile(Parameters &current_pars, std::vector<Pedigree> &ped1)
{

    std::string pedigreeFileName = current_pars.pedigree_file_path;
    std::ifstream inPedigreeFile;
    inPedigreeFile.open(pedigreeFileName, std::ios::in);

    if (!inPedigreeFile)
    {
        cerr << "Pedigree file could not be opened" << std::endl;
        exit(1);
    }

    std::string line;
    int id_count = 0;

    if (inPedigreeFile.is_open())
    {

        // Regex to match lines starting with one or more '*'
        std::regex pat_stars("^\\*+");

        bool foundDataStart = false;

        while (getline(inPedigreeFile, line)) {
            trim(line);
    
            // Skip all lines before the first line that starts with '**'
            if (!foundDataStart) {
                if (std::regex_match(line, pat_stars)) {
                    foundDataStart = true;
                }
                continue; // Skip everything before the marker line
            }
    
            // After marker: skip empty lines or comment lines
            if (line.empty()) {
                continue;
            }

                Pedigree temp1;
                id_count++;
                std::stringstream ss(line);
                ss >> temp1.subject_id_orig;
                ss >> temp1.father_id_orig;
                ss >> temp1.mother_id_orig;
                ss >> temp1.subject_sex;
                ss >> temp1.pheno;
                temp1.subject_id = to_string (id_count);
                ped1.push_back(temp1);
            
        }
    }

    inPedigreeFile.close();
}

void readDenseMarkerPositions(Parameters &current_pars, DenseMarkerPos &dense_marker_pos)
{

    std::string denseMarkerPosFileName = current_pars.dense_marker_pos_file_path;
    std::ifstream inDenseMarkerPosFile;
    inDenseMarkerPosFile.open(denseMarkerPosFileName, std::ios::in);

    if (!inDenseMarkerPosFile)
    {
        cerr << "Dense marker positions file could not be opened" << std::endl;
        exit(1);
    }

    std::string line;

    int cur_bp_pos;
    double cur_cM_pos;

    if (inDenseMarkerPosFile.is_open())
    {
        while (getline(inDenseMarkerPosFile, line))
        {
            std::stringstream ss(line);
            ss >> cur_bp_pos;
            ss >> cur_cM_pos;

            std::cout.setf(ios::fixed);
            std::cout.setf(ios::showpoint);
            std::cout.precision(6);

            dense_marker_pos.setVal(cur_bp_pos, cur_cM_pos);
        }
    }

    inDenseMarkerPosFile.close();
}

void readMeiosisFile(Parameters &current_pars, Pedigree_info &ped1,
MI_matrix *mi_matrix_ptr)
{

    std::string meiosisFileName = current_pars.meiosis_file_path;
    std::ifstream inMeiosisFile;
    inMeiosisFile.open(meiosisFileName, std::ios::in);

    if (!inMeiosisFile)
    {
        cerr << "Meiosis Indicators file could not be opened" << std::endl;
        exit(1);
    }

    int num_meiosis_lines = ped1.getNumMeiosisLines();

    int count_iter_num = 0;
    int count_m_lines = 0;
    std::string subject_id;
    int cur_meiosis_n = 0;
    int chr_origin = -1;
    int initial_chr = -1;
    int switches_num = -1;
    std::vector<int> cur_recomb_pos;
    int holder;
    string trash;

    std::string line;

    if (inMeiosisFile.is_open())
    {
        while (getline(inMeiosisFile, line))
        {
            ++count_m_lines;
            if (count_m_lines == 1)
            {
                ++count_iter_num;
            }

            std::string id_hold;
            std::stringstream ss(line); 
            ss >> id_hold;
            subject_id = ped1.convert2TempID(id_hold);
            ss >> trash;
            ss >> chr_origin;
            ss >> trash;
            ss >> initial_chr;
            ss >> switches_num;

            if (switches_num > 0)
            {
                while (ss)
                {
                    if (ss.rdbuf()->in_avail() == 0)
                        break;

                    ss >> holder;
                    cur_recomb_pos.push_back(holder);
                }
            }

            else

            {
                cur_recomb_pos.push_back(0);
            }

            int subject_index = ped1.getSubjectIndex(subject_id);
            cur_meiosis_n = ped1.getMIlineNum(subject_index, chr_origin);
            mi_matrix_ptr->setVal(count_iter_num, cur_meiosis_n, chr_origin,
            initial_chr, switches_num);

            mi_matrix_ptr->setRecombPos(count_iter_num, cur_meiosis_n, switches_num, cur_recomb_pos);
            
            cur_recomb_pos.clear();

            if (count_m_lines == num_meiosis_lines)
            {
                count_m_lines = 0;
            }
        }
    }

    inMeiosisFile.close();

}

int getNumSubWGS(Parameters &current_pars)
{

    std::string denseMarkerGenoTxtFileName = current_pars.dense_marker_geno_txt_file_path;
    std::ifstream inDenseMarkerGenoTxtFile;
    inDenseMarkerGenoTxtFile.open(denseMarkerGenoTxtFileName, std::ios::in);

    if (!inDenseMarkerGenoTxtFile)
    {
        cerr << "Dense marker genotypes file could not be opened" << std::endl;
        exit(1);
    }

    int num_ids_wgs = 0;
    std::string line;
    std::string temp1;
    int count2 = 0;

    if (inDenseMarkerGenoTxtFile.is_open())
    {
        if (getline(inDenseMarkerGenoTxtFile, line))
        {
            std::stringstream ss(line);

            while (ss >> temp1)
            {

                count2++;
            }
            goto next_step;
        }
    }

next_step:

    inDenseMarkerGenoTxtFile.close();

    num_ids_wgs = (count2 - 1) / 2;

    return num_ids_wgs;
}

void getChrNum(Parameters &current_pars)
{

    std::string denseMarkerGenoTxtFileName = current_pars.dense_marker_geno_txt_file_path;
    std::ifstream inDenseMarkerGenoTxtFile;
    inDenseMarkerGenoTxtFile.open(denseMarkerGenoTxtFileName, std::ios::in);

    if (!inDenseMarkerGenoTxtFile)
    {
        cerr << "Dense marker genotypes file could not be opened" << std::endl;
        exit(1);
    }

    std::string line;
    int count1 = 0;

    if (inDenseMarkerGenoTxtFile.is_open())
    {
        while (getline(inDenseMarkerGenoTxtFile, line))
        {
            count1++;

            if (count1 == 1)
            {

                continue;
            }

            else if (count1 == 2)

            {

                std::stringstream ss(line);
                std::string temp2;
                std::vector<std::string> temp3;
                std::vector<int> temp4;

                while (ss >> temp2)
                {
                    temp3.push_back(temp2);
                }

                temp4 = split(temp3[0], ':'); // temp4[0] = chr number
                current_pars.chr = temp4[0];

                goto next_step2;
            }
        }
    }

next_step2:

    inDenseMarkerGenoTxtFile.close();
}

void getInfoSubWGS(Parameters &current_pars, Pedigree_info &current_ped)
{

    std::string denseMarkerGenoTxtFileName = current_pars.dense_marker_geno_txt_file_path;
    std::ifstream inDenseMarkerGenoTxtFile;
    inDenseMarkerGenoTxtFile.open(denseMarkerGenoTxtFileName, std::ios::in);

    if (!inDenseMarkerGenoTxtFile)
    {
        cerr << "Dense marker genotypes file could not be opened" << std::endl;
        exit(1);
    }

    std::vector<int> col_index_conv; // provides subject index number for each column of data
    std::string line;
    std::string temp1;
    std::string junk;
    int count2 = 0;
    int pos_count = -1;

    if (inDenseMarkerGenoTxtFile.is_open())
    {
        if (getline(inDenseMarkerGenoTxtFile, line))
        {
            std::stringstream ss(line);

            ss >> junk;
            while (ss >> temp1)
            {

                count2++;

                std::string subject_id = current_ped.convert2TempID(temp1);
                int temp2 = current_ped.getSubjectIndex(subject_id);
                col_index_conv.push_back(temp2);

                if (count2 % 2 != 0)
                    continue;
                current_ped.cur_ped[temp2].wgs_data = 1;
                pos_count++;
                current_ped.cur_ped[temp2].mat_3D_pos = pos_count;
                pos_count++;
                current_ped.cur_ped[temp2].pat_3D_pos = pos_count;
            }
            goto next_step1;
        }
    }
next_step1:

    inDenseMarkerGenoTxtFile.close();
}

void readDenseMarkerGenotypes_txt(Parameters &current_pars, Pedigree_info &current_ped, DenseMarkerGenos &dense_genotypes, int start_line, int lines_to_read)
{

    std::string denseMarkerGenoTxtFileName = current_pars.dense_marker_geno_txt_file_path;
    std::ifstream inDenseMarkerGenoTxtFile;
    inDenseMarkerGenoTxtFile.open(denseMarkerGenoTxtFileName, std::ios::in);

    if (!inDenseMarkerGenoTxtFile)
    {
        cerr << "Dense marker genotypes file could not be opened" << std::endl;
        exit(1);
    }

    std::vector<int> col_index_conv; // provides subject index number for each column of data
    std::string line;
    std::string temp1;
    std::string junk;
    int count2 = 0;
    int pos_count = -1;

    if (inDenseMarkerGenoTxtFile.is_open())
    {
        if (getline(inDenseMarkerGenoTxtFile, line))
        {
            std::stringstream ss(line);

            ss >> junk;
            while (ss >> temp1)
            {

                count2++;

                std::string subject_id = current_ped.convert2TempID(temp1);
                int temp2 = current_ped.getSubjectIndex(subject_id);
                col_index_conv.push_back(temp2);

                if (count2 % 2 != 0)
                    continue;
                current_ped.cur_ped[temp2].wgs_data = 1;
                pos_count++;
                current_ped.cur_ped[temp2].mat_3D_pos = pos_count;
                pos_count++;
                current_ped.cur_ped[temp2].pat_3D_pos = pos_count;
            }
            goto next_step1;
        }
    }
next_step1:

    int end_line = (start_line + lines_to_read) - 1;

    int count1 = 0;

    if (inDenseMarkerGenoTxtFile.is_open())
    {
        while (getline(inDenseMarkerGenoTxtFile, line))
        {
            count1++;

            if (count1 < start_line)
            {
                continue;
            }

            else if (count1 >= start_line && count1 <= end_line)

            {

                std::stringstream ss(line);
                std::string temp2;
                std::vector<std::string> temp3;
                std::vector<int> temp4;

                while (ss >> temp2)
                {
                    temp3.push_back(temp2);
                }

                temp4 = split(temp3[0], ':'); // temp4[1] = bp position

                for (size_t i = 1; i < temp3.size(); i += 2)
                {
                    std::vector<std::string> temp5;
                    temp5.push_back(temp3[i]);
                    temp5.push_back(temp3[i + 1]);
                    dense_genotypes.setGeno(temp4[1], col_index_conv[i], temp5);
                }
            }

            else if (count1 > end_line)
            {
                goto next_step2;
            }
        }
    }

next_step2:

    inDenseMarkerGenoTxtFile.close();
}

void readFrameworkMarkerPositions(Parameters &current_pars, FrameworkMarkerPos &framework_marker_pos)
{

    std::string frameworkMarkerPosFileName = current_pars.framework_marker_pos_file_path;
    std::ifstream inFrameworkMarkerPosFile;
    inFrameworkMarkerPosFile.open(frameworkMarkerPosFileName, std::ios::in);

    if (!inFrameworkMarkerPosFile)
    {
        cerr << "Framework marker positions file could not be opened" << std::endl;
        exit(1);
    }

    std::string line;

    int count1 = 0;
    double cur_cM_pos;

    if (inFrameworkMarkerPosFile.is_open())
    {
        while (getline(inFrameworkMarkerPosFile, line))
        {
            count1++;
            std::stringstream ss(line);
            ss >> cur_cM_pos;

            std::cout.setf(ios::fixed);
            std::cout.setf(ios::showpoint);
            std::cout.precision(6);

            framework_marker_pos.setVal(cur_cM_pos);
        }
    }

    inFrameworkMarkerPosFile.close();
}

void readHaplotypeFile(Parameters &current_pars, int current_wind_1st_bp,
Pedigree_info *current_ped_ptr, HaploShare &wind_haplos)
{

    std::string haplotypeFileName = current_pars.haplotype_file_path;
    std::ifstream inHaplotypeFile;
    inHaplotypeFile.open(haplotypeFileName, std::ios::in);

    if (!inHaplotypeFile)
    {
        cerr << "File containing haplotypes could not be opened" << std::endl;
        exit(1);
    }

    Pedigree_info *ped1 = nullptr;
    ped1 = current_ped_ptr;

    int target_bp = current_wind_1st_bp;
    int current_bp = 0;
    int check1 = 0;
    std::string line;

    int count2 = 0;

    if (inHaplotypeFile.is_open())
    {
    LB1:;

        while (getline(inHaplotypeFile, line))
        {
            count2++;
            std::vector<char> st1(line.size() + 1);
            line.copy(st1.data(), line.size() + 1);
            st1[line.size()] = '\0';

            char *token1 = strtok(st1.data(), " ");

            while (token1 != NULL)
            {

                int count1 = 0;
                char *token2 = strtok(token1, ":-");
                while (token2 != NULL)
                {
                    count1++;

                    if (count1 == 2)
                    {

                        std::stringstream ss(token2);
                        ss >> current_bp;
                        break;
                    }

                    token2 = strtok(NULL, ":-");
                }

                token1 = strtok(NULL, " :");
                break;
            }

            if (current_bp == target_bp)
            {

                check1 = 7;

                std::vector<char> st2(line.size() + 1);
                line.copy(st2.data(), line.size() + 1);
                st2[line.size()] = '\0';

                char *token3 = strtok(st2.data(), " ");
                int count3 = 0;

                int chr_num_3D = -1;
                std::vector<int> haplo_seq;

                while (token3 != NULL)
                {
                    count3++;

                    if (count3 == 2)
                    {

                        string ind_chr = token3;
                        int cut_point = ind_chr.rfind("_");

                        int sub_id;
                        if (cut_point > -1)
                        {

                            std::string chr = ind_chr.substr(cut_point + 1, 1);

                            ind_chr.erase(cut_point);
                            std::string subject_id_tmp = ped1->convert2TempID(ind_chr);
                            sub_id = ped1->getSubjectIndex(subject_id_tmp);

                            if (chr == "0")
                            {
                                chr_num_3D = ped1->cur_ped[sub_id].mat_3D_pos;
                            }

                            if (chr == "1")
                            {
                                chr_num_3D = ped1->cur_ped[sub_id].pat_3D_pos;
                            }
                        }
                    }

                    if (count3 == 3)
                    {

                        std::string chr_haplo(token3);

                        for (size_t i = 0; i < chr_haplo.length(); ++i)
                        {

                            int haplo_al = (int)chr_haplo[i] - 48;

                            if (haplo_al == 3 || haplo_al == 7)
                            {

                                haplo_al = 0;
                            }

                            haplo_seq.push_back(haplo_al);
                        }
                    }

                    token3 = strtok(NULL, " ");
                }

                wind_haplos.setHaploData(chr_num_3D, haplo_seq);

                goto LB1;
            }

            if (check1 == 7)
            {
                goto LB2;
            }
        }
    }

LB2:;

    inHaplotypeFile.close();
}

void readCadidateGroup(Parameters &current_pars, Pedigree_info *current_ped_ptr)
{

    std::string candidateGrpFileName = current_pars.candidate_grp_file_path;
    std::ifstream inCandidateGrpFile;
    inCandidateGrpFile.open(candidateGrpFileName, std::ios::in);

    if (!inCandidateGrpFile)
    {
        cerr << "File containing IDs for subjects from the core group cannot not be opened" << std::endl;
        exit(1);
    }

    std::string line;
    std::vector<std::string> cadidate_grp;

    if (inCandidateGrpFile.is_open())
    {

        while (getline(inCandidateGrpFile, line))
        {

            std::stringstream ss(line);
            string temp1;

            while (ss >> temp1)
            {
                std::string subject_id = current_ped_ptr->convert2TempID(temp1);
                cadidate_grp.push_back(subject_id);
            }

            std::sort(cadidate_grp.begin(), cadidate_grp.end());
            current_pars.candidate_group = cadidate_grp;
        }
    }
}

void printCandidateGrp(ofstream *fgl_share_grp2, std::vector<int> &result1, Pedigree_info &current_ped)
{

    for (size_t i = 1; i < result1.size(); i++) // el[0] = inter #
    {
        std::string subject_id_orig = current_ped.convert2OrigID(std::to_string(result1[i]));
        *fgl_share_grp2 << subject_id_orig << " "; 
    }
    *fgl_share_grp2 << std::endl;
}

std::string getCurrentDir()
{

    char outDir[PATH_MAX];
    string path;
    if (getcwd(outDir, sizeof(outDir)) != NULL)
    {
        return path = outDir;
    }
    else
    {
        perror("getcwd() error");
    }

    return ("NA");
}

ofstream *createOutStream(std::string file_path)
{

    ofstream *fileStream = new ofstream;
    fileStream->open(file_path);

    if (!(*fileStream))
    {
        cerr << "Cannot open" << file_path << "output file for writing data in." << std::endl;
        exit(1);
    }
    return fileStream;
}


void show_help() 
{

    std::string message = R"(
usage: HaploGI [option] 

            options:
                        --version 
                        displays program version number

                        --help 
                        displays program run options

usage: HaploGI [run_option] parameter_file

        run options:
                        --haplotyping  
                        run option to perform haplotyping
                        of whole genome sequence (WGS) data

                        --haplosharing 
                        run option to determine haplotype sharing among
                        subjects with whole genome sequence (WGS) data

                        --full 
                        run option to perform haplotyping followed by
                        by determination of haplotype sharing
                        
     parameter_file: 
                        path to a parameter file
    )";
        
    std::cout << message << "\n" << std::endl;
        
}

void show_program_info() 
{

std::string message = R"(
HaploGI - version: 1.0.0
© Rafael Nafikov 2025

Licenced under GNU General Public License (GPLv3), 
see <https://www.gnu.org/licenses/>

HaploGI: Haplotyping Given Inheritance

Comments and suggestions should be emailed to:
	nrafscience@gmail.com
)";

std::cout << message << "\n" << std::endl;
   
}

class TeeBuf : public std::streambuf {
    public:
        TeeBuf(std::streambuf* buf1, std::streambuf* buf2)
            : m_buf1(buf1), m_buf2(buf2) {}
    
    protected:
        virtual int overflow(int c) override {
            if (c == EOF) return !EOF;
            int const r1 = m_buf1->sputc(c);
            int const r2 = m_buf2->sputc(c);
            return (r1 == EOF || r2 == EOF) ? EOF : c;
        }
    
        virtual int sync() override {
            int const r1 = m_buf1->pubsync();
            int const r2 = m_buf2->pubsync();
            return (r1 == 0 && r2 == 0) ? 0 : -1;
        }
    
    private:
        std::streambuf* m_buf1;
        std::streambuf* m_buf2;
    };



///  MAIN 
int main(int argc, char *argv[])
{


/// determining run type
    int run_type = 0;

    std::string parameterFileName;

    if (argc == 1)
    {
        show_program_info();
    }

    for (int i = 1; i < argc; ++i)
    {

        if (std::string(argv[i]) == "--haplotyping")
        {
            if (argc < 3)
            {
                cerr << "Run option and/or parameter file are not specified correctly!" << std::endl;
                exit(1);
            }

            else if (i + 1 < argc)
            {
                parameterFileName = argv[++i];
       
                run_type = 1;
            }
        }

        else if (std::string(argv[i]) == "--haplosharing")
        {
            if (argc < 3)
            {
                cerr << "Run option and/or parameter file are not specified correctly!" << std::endl;
                exit(1);
            }
            
            else if (i + 1 < argc)
            {
                parameterFileName = argv[++i];

                run_type = 2;
            }
        }

        else if (std::string(argv[i]) == "--full")
        {
            if (argc < 3)
            {
                cerr << "Run option and/or parameter file are not specified correctly!" << std::endl;
                exit(1);
            }

            else if (i + 1 < argc)
            {
                parameterFileName = argv[++i];

                run_type = 3;
            }
        }

        else if (std::string(argv[i]) == "--version" || std::string(argv[i]) == "-v" )
        {
            show_program_info();
        }

        else if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h" )
        {
            show_program_info();
            show_help();
        }

        else
        {
            cerr << "Run option and/or parameter file are not specified correctly!" << std::endl;
            exit(1);
        }
    }

    /// HAPLOTYPING ONLY RUN
    if (run_type == 1)
    {
        
        auto run_start = std::chrono::high_resolution_clock::now();

        /// reading par file
        Parameters current_pars;
        current_pars.run_type = 1;
        readParameterFile(parameterFileName, current_pars);
        Parameters *current_pars_ptr = &current_pars;
        int window_size = current_pars.window_size;
        std::string out_dir = current_pars.output_dir_path;


        /// logging std::cout and cerr messages
        std::string file_name4 = out_dir + "/log.txt";
        std::ofstream logStream (file_name4);
       
        std::streambuf* coutBuf = std::cout.rdbuf();
        std::streambuf* cerrBuf = std::cerr.rdbuf();
 
        TeeBuf teeBufOut(coutBuf, logStream.rdbuf());
        TeeBuf teeBufErr(cerrBuf, logStream.rdbuf());
        
        std::cout.rdbuf(&teeBufOut);
        std::cerr.rdbuf(&teeBufErr);
 
        show_program_info();
        std::cout << "Begining of haplotyping run" << "\n" << std::endl;

        /// current date and time
        auto now = std::chrono::system_clock::now();
        std::time_t now_time = std::chrono::system_clock::to_time_t(now);
        std::cout << "Current date and time: " << std::ctime(&now_time) << std::endl;


        /// getting seed # from par file and passing it to rand number generator
        int seed = current_pars_ptr->seed;
        std::mt19937 rng(seed);// always pass by reference
        std::cout << "Seed number used is: " << seed << std::endl;
        
        /// printout of additional run conditions
        std::cout << "Number of iterations in meiosis indicators (MI) file is: " << current_pars.num_iter << std::endl;
        std::cout << "Boundaries of linkage region are: " << current_pars.linkage_region[0] << " " << current_pars.linkage_region[1] << std::endl;
        std::cout << "Linkage marker with maxLOD score is at : " << current_pars.maxLODmarker << " cM" << "\n" << std::endl;

       
        /// generating filestreams and output files
        std::string file_name1 = out_dir + "/core_set_of_cases.txt";
        std::ofstream *fglShareStream2 = createOutStream(file_name1);
        
        std::string file_name2 = out_dir + "/allele2FGLinconsistency.txt";
        std::ofstream *consistencyStream = createOutStream(file_name2);
        *consistencyStream << "bp_position" << std::endl;
        
        std::string file_name3 = out_dir + "/haplotype_sequences.txt";
        std::ofstream *haploStream2 = createOutStream(file_name3);
       
        
        /// getting chr# from WGS data file
        getChrNum(current_pars);
        int chrNumIn = current_pars.chr;

      
        /// reading and setting pedigree file info
        std::cout << "Reading pedigree file" << "\n" << std::endl;
        std::vector<Pedigree> ped1;
        readPedigreeFile(current_pars, ped1);
        Pedigree_info current_ped(ped1);
        Pedigree_info *current_ped_ptr = &current_ped;
        current_ped.generateTempIDs();
        current_ped.genMI_FGLrelInfo();
        

        /// reading MI file
        std::cout << "Reading MI file" << "\n" << std::endl;
        MI_matrix mi_matrix(current_pars.num_iter, current_ped.getNumMeiosisLines());
        MI_matrix *mi_matrix_ptr = &mi_matrix;
        readMeiosisFile(current_pars, current_ped, mi_matrix_ptr);


        /// reading genomic positions for SNPs (linkage markers)
        std::cout << "Reading genomic positions for linkage markers" << "\n" << std::endl;
        FrameworkMarkerPos framework_marker_pos;
        FrameworkMarkerPos *framework_marker_pos_ptr = &framework_marker_pos;
        readFrameworkMarkerPositions(current_pars, framework_marker_pos);


        /// reading SNV positions file
        std::cout << "Reading genomic positions for SNVs" << "\n" << std::endl;
        DenseMarkerPos dense_markers_pos(window_size);
        readDenseMarkerPositions(current_pars, dense_markers_pos);
        dense_markers_pos.makeWindows();
        DenseMarkerPos *dense_markers_pos_ptr = &dense_markers_pos;


        /// creates conversion system to read genotypes of subjects with WGS data from 3D matrix
        getInfoSubWGS(current_pars, current_ped);


        /// determination of FGL sharing with identification of the core set of cases and iteration # for phasing
        std::cout << "Identifying core set of cases and iteration # for phasing" << "\n" << std::endl;
        FGLshareGroup2 fgl_share_group2(current_ped_ptr, current_pars_ptr, mi_matrix_ptr,
        framework_marker_pos_ptr, dense_markers_pos_ptr, rng);
        std::vector<int> resultShareGrp = fgl_share_group2.evalFGLshareMain();
        printCandidateGrp(fglShareStream2, resultShareGrp, current_ped);


        /// determining total # of genomic windows based on # of SNVs in WGS data
        int t_num_windows = dense_markers_pos.getNumWindows();

        std::cout << "Starting phasing and haplotyping of WGS data" << "\n" << std::endl;
        std::cout << "Total number of genomic windows to be used is: " << t_num_windows << "\n" << std::endl; 
        std::cout << "Reporting progress every 100 genomic windows" << "\n" << std::endl;

        /// START OF GENOMIC WINDOWS LOOP
        for (int cur_wind_num = 1; cur_wind_num <= t_num_windows; cur_wind_num++)
        {

            if (cur_wind_num % 100 == 0) 
            {
                std::cout << "Completed processing " << cur_wind_num << " genomic windows" << std::endl;
            } 

            /// getting WGS data for a genomic window
            int start_line = (window_size * (cur_wind_num - 1)) + 1;
            DenseMarkerGenos dense_genotypes;
            readDenseMarkerGenotypes_txt(current_pars, current_ped, dense_genotypes, start_line, window_size);
            DenseMarkerGenos *dense_genotypes_ptr = &dense_genotypes;


            /// getting bp positions for SNVs in a genomic window
            std::vector<int> cur_wind_bps;
            cur_wind_bps = dense_markers_pos.getWindowBps(cur_wind_num);
            std::vector<int> *cur_wind_bps_ptr = &cur_wind_bps;
            

            /// creating data structure to store phased haplotypes for specific window
            MatrixHaplo haplo_window2(window_size, cur_wind_bps_ptr, current_ped_ptr, chrNumIn);
            MatrixHaplo *haplo_window2_ptr = &haplo_window2;

           
            /// START OF BP POSITIONS LOOP within a genomic window
            for (std::vector<int>::iterator iter1 = cur_wind_bps.begin();
            iter1 != cur_wind_bps.end(); iter1++)
            {
                
                /// FGL matrix for specific iteration and bp position within window
                FGL_matrix fgl_matrix(current_pars.num_iter, current_ped.getNumSubjects(),
                current_ped_ptr, mi_matrix_ptr, framework_marker_pos_ptr, rng, dense_markers_pos_ptr);
                fgl_matrix.setFGL(*iter1);
                FGL_matrix *fgl_matrix_ptr = &fgl_matrix;

                /// determines sequential var number within window (max = window size)
                int var_pos = std::distance(cur_wind_bps.begin(), iter1); // starts from 0

                /// passing iter number for phasing
                int iter_num2phase = resultShareGrp[0]; // passes specific iter# starting from 1 to total # of iter

                Phase2 phased_genos2(fgl_matrix_ptr, dense_genotypes_ptr, current_ped_ptr,
                haplo_window2_ptr, *iter1, var_pos, iter_num2phase, consistencyStream);

                phased_genos2.phaseGenos();

                
             

            } /// END OF BP POSITIONS LOOP within a genomic window

           
            /// outputting phased haplotype sequences for specific window into the file
            haplo_window2.printOutResult(haploStream2);


        } /// END of GENOMIC WINDOWS LOOP


        fglShareStream2->close();
        consistencyStream->close();
        haploStream2->close();
        
        // MI 0 = Maternal, 1 = Paternal

        std::cout << "\n" << "Haplotyping run is complete!" << "\n" << std::endl;
        auto run_end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(run_end - run_start).count();
        // Convert to hours, minutes, seconds
        int hours = duration / 3600;
        int minutes = (duration % 3600) / 60;
        int seconds = duration % 60;
        std::cout << "Program runtime: " << hours << "h " << minutes << "m " << seconds << "s" << "\n" << std::endl;



        logStream.close();

        // restoring buffers
        std::cout.rdbuf(coutBuf);
        std::cerr.rdbuf(cerrBuf);
    } /// END of run_type 1


    /// HAPLOSHARING ONLY RUN
    else if (run_type == 2)
    {
        
        auto run_start = std::chrono::high_resolution_clock::now();
        
        /// reading par file
        Parameters current_pars;
        current_pars.run_type = 2;
        readParameterFile(parameterFileName, current_pars);
        Parameters *current_pars_ptr = &current_pars;
        int window_size = current_pars.window_size;
        std::string out_dir = current_pars.output_dir_path;


        /// logging std::cout and cerr messages
        std::string file_name4 = out_dir + "/log.txt";
        std::ofstream logStream (file_name4);
        
        std::streambuf* coutBuf = std::cout.rdbuf();
        std::streambuf* cerrBuf = std::cerr.rdbuf();
  
        TeeBuf teeBufOut(coutBuf, logStream.rdbuf());
        TeeBuf teeBufErr(cerrBuf, logStream.rdbuf());
         
        std::cout.rdbuf(&teeBufOut);
        std::cerr.rdbuf(&teeBufErr);

        show_program_info();
        std::cout << "Begining of haplosharing run" << "\n" << std::endl;

        /// current date and time
        auto now = std::chrono::system_clock::now();
        std::time_t now_time = std::chrono::system_clock::to_time_t(now);
        std::cout << "Current date and time: " << std::ctime(&now_time) << std::endl;


        /// getting seend # from par file and passing it to rand number generator
        int seed = current_pars_ptr->seed;
        std::mt19937 rng(seed);// always pass by reference
        std::cout << "Seed number used is: " << seed << std::endl;
       
        /// printout of additional run conditions
        std::cout << "Boundaries of linkage region are: " << current_pars.linkage_region[0] << " " << current_pars.linkage_region[1] << std::endl;
        std::cout << "Linkage marker with maxLOD score is at : " << current_pars.maxLODmarker << " cM" << "\n" << std::endl;

        /// generating filestreams and output files
        std::string file_name1 = out_dir + "/shared_haplotype_temp.txt";
        std::ofstream *sharedHaploStream = createOutStream(file_name1);
        
        std::string file_name2 = out_dir + "/haplotype_sharing_patterns.txt";
        std::ofstream *haploSharingStream = createOutStream(file_name2);
        
        std::string file_name3 = out_dir + "/shared_haplotype.txt";
        std::ofstream *sharedHaploStreamFinal = createOutStream(file_name3);

        /// getting chr# from WGS data file
        getChrNum(current_pars);


        /// reading and setting pedigree file info
        std::cout << "Reading pedigree file" << "\n" << std::endl;
        std::vector<Pedigree> ped1;
        readPedigreeFile(current_pars, ped1);
        Pedigree_info current_ped(ped1);
        Pedigree_info *current_ped_ptr = &current_ped;
        current_ped.generateTempIDs();

        /// reading supplied the core set of cases from file
        std::cout << "Reading IDs of cases from core set" << "\n" << std::endl;
        readCadidateGroup(current_pars, current_ped_ptr);


        /// reading SNV positions file
        std::cout << "Reading genomic positions for SNVs" << "\n" << std::endl;
        DenseMarkerPos dense_markers_pos(window_size);
        readDenseMarkerPositions(current_pars, dense_markers_pos);
        dense_markers_pos.makeWindows();
        DenseMarkerPos *dense_markers_pos_ptr = &dense_markers_pos;


        /// determining total # of genomic windows based on # of SNVs in WGS data
        int t_num_windows = dense_markers_pos.getNumWindows();


        /// creates conversion system to read genotypes of subjects with WGS data from 3D matrix
        int numSubWGS = getNumSubWGS(current_pars);
        
        
        /// creates data structure to store results of haplotype sharing determination in ROI
        IntegrationHpSh global_hp_sharing(numSubWGS, t_num_windows, current_ped_ptr, haploSharingStream, sharedHaploStreamFinal, current_pars_ptr, dense_markers_pos_ptr);
        IntegrationHpSh *global_hp_sharing_ptr = &global_hp_sharing;

        std::cout << "Starting determination of haplotype sharing" << "\n" << std::endl;
        std::cout << "Total number of genomic windows to be used is: " << t_num_windows << "\n" << std::endl; 
        std::cout << "Reporting progress every 100 genomic windows" << "\n" << std::endl;

        /// START OF GENOMIC WINDOWS LOOP
        for (int cur_wind_num = 1; cur_wind_num <= t_num_windows; cur_wind_num++)
        {
             
            if (cur_wind_num % 100 == 0) 
            {
                std::cout << "Completed processing " << cur_wind_num << " genomic windows" << std::endl;
            } 


            /// getting WGS data for a genomic window
            int start_line = (window_size * (cur_wind_num - 1)) + 1;
            DenseMarkerGenos dense_genotypes;
            readDenseMarkerGenotypes_txt(current_pars, current_ped, dense_genotypes, start_line, window_size);
            DenseMarkerGenos *dense_genotypes_ptr = &dense_genotypes;


            /// getting bp positions for SNVs in a genomic window
            std::vector<int> cur_wind_bps;
            cur_wind_bps = dense_markers_pos.getWindowBps(cur_wind_num);
            std::vector<int> *cur_wind_bps_ptr = &cur_wind_bps;


            /// determination of haplotype sharing within a genomic window
            HaploShare wind_haplos(window_size, cur_wind_num, current_ped_ptr,
            global_hp_sharing_ptr, current_pars, dense_genotypes_ptr, cur_wind_bps_ptr, sharedHaploStream);
            readHaplotypeFile(current_pars, cur_wind_bps[0], current_ped_ptr, wind_haplos);
            wind_haplos.analyzeHpShP1();

          
        } /// END OF GENOMIC WINDOWS LOOP


        /// summarizing haplosharing results across ROI
        int sharedHpFileFate = global_hp_sharing.integrateHpSh();


        /// start of handaling of temporary file containing shared haplotype sequence
        int n1 = file_name1.length();
        int n3 = file_name3.length();
        std::vector<char> file1(n1 + 1);
        std::vector<char> file3(n3 + 1);
        std::strcpy(file1.data(), file_name1.c_str());
        std::strcpy(file3.data(), file_name3.c_str());

        if (sharedHpFileFate == -1)
        {

            if (remove(file3.data()) != 0)
            {
                std::cout << "shared_haplotype.txt cannot be deleted! Please delete this file manually";
            }

            if (rename(file1.data(), file3.data()) != 0)
            {
                std::cout << "shared_haplotype_temp.txt cannot be renamed into shared_haplotype.txt! Please rename this file manually!";
            }
        }

        else

        {

            if (remove(file1.data()) != 0)
            {
                std::cout << "shared_haplotype_temp.txt cannot be deleted! Please disregard this file!";
            }
        }
        /// end of handaling of temporary file

        sharedHaploStream->close();
        haploSharingStream->close();
        sharedHaploStreamFinal->close();

        std::cout << "\n" << "Haplosharing run is complete!" << "\n" << std::endl;
        auto run_end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(run_end - run_start).count();
        // Convert to hours, minutes, seconds
        int hours = duration / 3600;
        int minutes = (duration % 3600) / 60;
        int seconds = duration % 60;
        std::cout << "Program runtime: " << hours << "h " << minutes << "m " << seconds << "s" << "\n" << std::endl;

        logStream.close();

        // restoring buffers
        std::cout.rdbuf(coutBuf);
        std::cerr.rdbuf(cerrBuf);
    } /// END of RUN_TYPE 2

   

    /// HAPLOTYPING AND HAPLOSHARING RUN
    else if (run_type == 3)
    {

        auto run_start = std::chrono::high_resolution_clock::now();

        /// reading par file
        Parameters current_pars;
        current_pars.run_type = 3;
        readParameterFile(parameterFileName, current_pars);
        Parameters *current_pars_ptr = &current_pars;
        int window_size = current_pars.window_size;
        std::string out_dir = current_pars.output_dir_path;


        /// logging std::cout and cerr messages
        std::string file_name7 = out_dir + "/log.txt";
        std::ofstream logStream (file_name7);
        
        std::streambuf* coutBuf = std::cout.rdbuf();
        std::streambuf* cerrBuf = std::cerr.rdbuf();
  
        TeeBuf teeBufOut(coutBuf, logStream.rdbuf());
        TeeBuf teeBufErr(cerrBuf, logStream.rdbuf());
         
        std::cout.rdbuf(&teeBufOut);
        std::cerr.rdbuf(&teeBufErr);

        show_program_info();
        std::cout << "Begining of full run" << "\n" << std::endl;

        /// current date and time
        auto now = std::chrono::system_clock::now();
        std::time_t now_time = std::chrono::system_clock::to_time_t(now);
        std::cout << "Current date and time: " << std::ctime(&now_time) << std::endl;


        /// getting seend # from par file and passing it to rand number generator
        int seed = current_pars_ptr->seed;
        std::mt19937 rng(seed);// always pass by reference
        std::cout << "Seed number used is: " << seed << std::endl;

        /// printout of additional run conditions
        std::cout << "Number of iterations in meiosis indicators (MI) file is: " << current_pars.num_iter << std::endl;
        std::cout << "Boundaries of linkage region are: " << current_pars.linkage_region[0] << " " << current_pars.linkage_region[1] << std::endl;
        std::cout << "Linkage marker with maxLOD score is at : " << current_pars.maxLODmarker << " cM" << "\n" << std::endl;

        /// generating filestreams and output files
        std::string file_name1 = out_dir + "/core_set_of_cases.txt";
        std::ofstream *fglShareStream2 = createOutStream(file_name1);

        std::string file_name2 = out_dir + "/allele2FGLinconsistency.txt";
        std::ofstream *consistencyStream = createOutStream(file_name2);
        *consistencyStream << "bp_position" << std::endl;

        std::string file_name3 = out_dir + "/haplotype_sequences.txt";
        std::ofstream *haploStream2 = createOutStream(file_name3);
        
        
        /// getting chr# from WGS data file
        getChrNum(current_pars);
        int chrNumIn = current_pars.chr;


        /// reading and setting pedigree file info
        std::cout << "Reading pedigree file" << "\n" << std::endl;
        std::vector<Pedigree> ped1;
        readPedigreeFile(current_pars, ped1);
        Pedigree_info current_ped(ped1);
        Pedigree_info *current_ped_ptr = &current_ped;
        current_ped.generateTempIDs();
        current_ped.genMI_FGLrelInfo();


        /// reading MI file
        std::cout << "Reading MI file" << "\n" << std::endl;
        MI_matrix mi_matrix(current_pars.num_iter, current_ped.getNumMeiosisLines());
        MI_matrix *mi_matrix_ptr = &mi_matrix;
        readMeiosisFile(current_pars, current_ped, mi_matrix_ptr);


        /// reading genomic positions for SNPs (linkage markers)
        std::cout << "Reading genomic positions for linkage markers" << "\n" << std::endl;
        FrameworkMarkerPos framework_marker_pos;
        FrameworkMarkerPos *framework_marker_pos_ptr = &framework_marker_pos;
        readFrameworkMarkerPositions(current_pars, framework_marker_pos);


        /// reading SNV positions file
        std::cout << "Reading genomic positions for SNVs" << "\n" << std::endl;
        DenseMarkerPos dense_markers_pos(window_size);
        readDenseMarkerPositions(current_pars, dense_markers_pos);
        dense_markers_pos.makeWindows();
        DenseMarkerPos *dense_markers_pos_ptr = &dense_markers_pos;


        /// creates conversion system to read genotypes of subjects with WGS data from 3D matrix
        getInfoSubWGS(current_pars, current_ped);


        /// determination of FGL sharing with identification of the core set of cases and iteration # for phasing
        std::cout << "Identifying core set of cases and iteration # for phasing" << "\n" << std::endl;
        FGLshareGroup2 fgl_share_group2(current_ped_ptr, current_pars_ptr, mi_matrix_ptr,
        framework_marker_pos_ptr, dense_markers_pos_ptr, rng);
        std::vector<int> resultShareGrp = fgl_share_group2.evalFGLshareMain();
        printCandidateGrp(fglShareStream2, resultShareGrp, current_ped);


        /// determining total # of genomic windows based on # of SNVs in WGS data
        int t_num_windows = dense_markers_pos.getNumWindows();
        
        std::cout << "Starting phasing and haplotyping of WGS data" << "\n" << std::endl;
        std::cout << "Total number of genomic windows to be used is: " << t_num_windows << "\n" << std::endl; 
        std::cout << "Reporting progress every 100 genomic windows" << "\n" << std::endl;

        /// START OF GENOMIC WINDOWS LOOP
        for (int cur_wind_num = 1; cur_wind_num <= t_num_windows; cur_wind_num++)
        {

            if (cur_wind_num % 100 == 0) 
            {
                std::cout << "Completed processing " << cur_wind_num << " genomic windows" << std::endl;
            } 

            
            /// getting WGS data for a genomic window
            int start_line = (window_size * (cur_wind_num - 1)) + 1;
            DenseMarkerGenos dense_genotypes;
            readDenseMarkerGenotypes_txt(current_pars, current_ped, dense_genotypes, start_line, window_size);
            DenseMarkerGenos *dense_genotypes_ptr = &dense_genotypes;


            /// getting bp positions for SNVs in a genomic window
            std::vector<int> cur_wind_bps;
            cur_wind_bps = dense_markers_pos.getWindowBps(cur_wind_num);
            std::vector<int> *cur_wind_bps_ptr = &cur_wind_bps;
            
            /// creating data structure to store phased haplotypes for specific window
            MatrixHaplo haplo_window2(window_size, cur_wind_bps_ptr, current_ped_ptr, chrNumIn);
            MatrixHaplo *haplo_window2_ptr = &haplo_window2;


            /// START OF BP POSITIONS LOOP within a genomic window
            for (std::vector<int>::iterator iter1 = cur_wind_bps.begin();
            iter1 != cur_wind_bps.end(); iter1++)
            {


                /// FGL matrix for specific iteration and bp position within window
                FGL_matrix fgl_matrix(current_pars.num_iter, current_ped.getNumSubjects(),
                current_ped_ptr, mi_matrix_ptr, framework_marker_pos_ptr, rng, dense_markers_pos_ptr);
                fgl_matrix.setFGL(*iter1);
                FGL_matrix *fgl_matrix_ptr = &fgl_matrix;


                /// determines sequential var number within window (max = window size)
                int var_pos = std::distance(cur_wind_bps.begin(), iter1); // starts from 0

                /// passing iter number for phasing
                int iter_num2phase = resultShareGrp[0]; // passes specific iter# starting from 1 to total # of iter

                Phase2 phased_genos2(fgl_matrix_ptr, dense_genotypes_ptr, current_ped_ptr,
                haplo_window2_ptr, *iter1, var_pos, iter_num2phase, consistencyStream);

                phased_genos2.phaseGenos();


    
            } /// END OF BP POSITIONS LOOP within a genomic window

            /// outputting phased haplotype sequences for a genomic window into the file
            haplo_window2.printOutResult(haploStream2);



        } /// END OF GENOMIC WINDOWS LOOP
         
        fglShareStream2->close();
        consistencyStream->close();// not used now
        haploStream2->close();
      

        /// MI 0 = Maternal, 1 = Paternal
    
        /// END OF PART 1 - HAPLOTYPING
        ///////////////////////////////////////////////////////
        
        /// PART 2 - HAPLOSHARING


        /// generating filestreams and output files
        std::string file_name4 = out_dir + "/shared_haplotype_temp.txt";
        std::ofstream *sharedHaploStream = createOutStream(file_name4);
        
        std::string file_name5 = out_dir + "/haplotype_sharing_patterns.txt";
        std::ofstream *haploSharingStream = createOutStream(file_name5);
        
        std::string file_name6 = out_dir + "/shared_haplotype.txt";
        std::ofstream *sharedHaploStreamFinal = createOutStream(file_name6);

        /// reading, determined in part 1, the core set of cases from file
        readCadidateGroup(current_pars, current_ped_ptr);


        /// returns # of subjects with WGS data
        int numSubWGS = getNumSubWGS(current_pars);

        /// creates data structure to store results of haplotype sharing determination in ROI
        IntegrationHpSh global_hp_sharing(numSubWGS, t_num_windows, current_ped_ptr, haploSharingStream, sharedHaploStreamFinal, current_pars_ptr, dense_markers_pos_ptr);
        IntegrationHpSh *global_hp_sharing_ptr = &global_hp_sharing;

        std::cout << "\n" << "Starting determination of haplotype sharing" << "\n" << std::endl;
        std::cout << "Total number of genomic windows to be used is: " << t_num_windows << "\n" << std::endl; 
        std::cout << "Reporting progress every 100 genomic windows" << "\n" << std::endl;

        /// START OF GENOMIC WINDOWS LOOP
        for (int cur_wind_num = 1; cur_wind_num <= t_num_windows; cur_wind_num++)
        {

            if (cur_wind_num % 100 == 0) 
            {
                std::cout << "Completed processing " << cur_wind_num << " genomic windows" << std::endl;
            } 

            /// getting WGS data for a genomic window
            int start_line = (window_size * (cur_wind_num - 1)) + 1;
            DenseMarkerGenos dense_genotypes;
            readDenseMarkerGenotypes_txt(current_pars, current_ped, dense_genotypes, start_line, window_size);
            DenseMarkerGenos *dense_genotypes_ptr = &dense_genotypes;


            /// getting bp positions for SNVs in a genomic window
            std::vector<int> cur_wind_bps;
            cur_wind_bps = dense_markers_pos.getWindowBps(cur_wind_num);
            std::vector<int> *cur_wind_bps_ptr = &cur_wind_bps;


            /// determination of haplotype sharing within a genomic window
            HaploShare wind_haplos(window_size, cur_wind_num, current_ped_ptr,
            global_hp_sharing_ptr, current_pars, dense_genotypes_ptr, cur_wind_bps_ptr, sharedHaploStream);
            readHaplotypeFile(current_pars, cur_wind_bps[0], current_ped_ptr, wind_haplos);
            wind_haplos.analyzeHpShP1();


        } /// END OF GENOMIC WINDOWS LOOP



        /// summarizing haplosharing results across ROI
        int sharedHpFileFate = global_hp_sharing.integrateHpSh();


        /// start of handaling of temporary file containing shared haplotype sequence
        int n1 = file_name4.length();
        int n3 = file_name6.length();
        std::vector<char> file1(n1 + 1);
        std::vector<char> file3(n3 + 1);

        std::strcpy(file1.data(), file_name4.c_str());
        std::strcpy(file3.data(), file_name6.c_str());

        if (sharedHpFileFate == -1)
        {

            if (remove(file3.data()) != 0)
            {
                std::cout << "shared_haplotype.txt cannot be deleted! Please delete this file manually";
            }

            if (rename(file1.data(), file3.data()) != 0)
            {
                std::cout << "shared_haplotype_temp.txt cannot be renamed into shared_haplotype.txt! Please rename this file manually!";
            }
        }

        else

        {

            if (remove(file1.data()) != 0)
            {
                std::cout << "shared_haplotype_temp.txt cannot be deleted! Please disregard this file!";
            }
        }
        /// end of handaling of temporary file

        sharedHaploStream->close();
        haploSharingStream->close();
        sharedHaploStreamFinal->close();

        std::cout << "\n" << "Full run is complete!" << "\n" << std::endl;
        auto run_end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(run_end - run_start).count();
        // Convert to hours, minutes, seconds
        int hours = duration / 3600;
        int minutes = (duration % 3600) / 60;
        int seconds = duration % 60;
        std::cout << "Program runtime: " << hours << "h " << minutes << "m " << seconds << "s" << "\n" << std::endl;

        logStream.close();

        // restoring buffers
        std::cout.rdbuf(coutBuf);
        std::cerr.rdbuf(cerrBuf);
    } /// END of RUN_TYPE 3

    return (0);
} /// END OF MAIN
