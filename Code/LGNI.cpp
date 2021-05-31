#include <iostream> 
#include <list> 
#include <algorithm>
#include <fstream>
#include <map>
#include <queue>
#include <time.h>
#include <cstring>
#include <vector>
#include <valarray>
#include <numeric> 
#include <algorithm> 
using namespace std; 
#include <chrono> 
#include <pthread.h>
#include <omp.h>
#include <stdlib.h>
using namespace std::chrono; 
  
template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
  // reverse(idx.begin(), idx.end());
  return idx;
}

ofstream ff;

class C_Node
{
public:
    bool is_controlled;
    float p_infect;
    float cost;
    C_Node()
    {
        is_controlled = false;
        p_infect = 0;
        cost = 0;
    }
};
class C_Edge
{
public:
    int v;
    float w;
    C_Edge(int v, float w)
    {
        this->v = v;
        this->w = w;
    }
};
class C_Visit
{
public:
    int u;
    int k;
    int v;
    C_Visit(int u, int k, int v)
    {
        this->u = u;
        this->k = k;
        this->v = v;
    }
};

class C_Ret
{
public:
    float h;
    float *g;
    float *hh;
    float *hhh;
    C_Ret(int V)
    {
        h = 0;
        g = new float[V];
        hh = new float[V];
        hhh = new float[V];
        // memset(g, 1, sizeof(g));
        // memset(hh, 0, sizeof(hh));
        for (int i = 0; i < V; i++)
        {
            this->g[i] = 1;
            this->hh[i] = 0;
            this->hhh[i] = 0;
        }
    }
    C_Ret()
    {
        g = new float[2];
        hh = new float[2];
        hhh = new float[2];
    }
    ~C_Ret()
    {
        delete[] this->g;
        delete[] this->hh;
        delete[] this->hhh;
    }
};

class Graph 
{
public:  
    int V;    
    int num_affected;
    int num_unaffected;
    float p_to_l;
    float l_to_p;
    list<C_Edge> *adj; 
    list<C_Edge> *re_adj;
    C_Node *nodes;
  
    Graph(int V, int num_affected, int num_unaffected, float p_to_l, float l_to_p);   // Constructor 
    ~Graph(); 
    void addEdge(int u, int v, float w);
    void addNode(int u, C_Node c);
    list<int> run_select(list<int> count, int k, float remaining_cost, float &minus_value, int step);
    C_Ret iterative_calculate(int u, int k, bool is_infected, float p_to_l, float l_to_p);
    float dfs(int cur_node, int step, const int uaf_count[], float** values, bool** flags);
}; 
  
// Method to print connected components in an 
// undirected graph 

Graph::Graph(int V, int num_affected, int num_unaffected, float p_to_l, float l_to_p)
{ 
    this->V = V; 
    this->num_affected = num_affected;
    this->num_unaffected = num_unaffected;
    this->p_to_l = p_to_l;
    this->l_to_p = l_to_p;
    this->adj = new list<C_Edge>[V]; 
    this->re_adj = new list<C_Edge>[V];
    this->nodes = new C_Node[V];
} 

Graph::~Graph() 
{ 
    delete[] adj; 
    delete[] re_adj;
} 
  
void Graph::addEdge(int u, int v, float w) 
{
    C_Edge c(v, w); 
    adj[u].push_back(c); 
    C_Edge d(u, w);
    re_adj[v].push_back(d);
} 

void Graph::addNode(int u, C_Node c)
{
    nodes[u] = c;
}

float Graph::dfs(int cur_node, int step, const int uaf_count[], float** values, bool** flags)
{
    float result = uaf_count[cur_node];
    if (step == 0)
        return result;
    for (list<C_Edge>::iterator j = this->re_adj[cur_node].begin(); j != this->re_adj[cur_node].end(); ++j)
    {
        if (!flags[j->v][step-1])
        {
            values[j->v][step-1] = dfs(j->v, step-1, uaf_count, values, flags);
            flags[j->v][step-1] = true;
        }
        result += j->w * values[j->v][step-1];
    }
    return result;
}

list<int> Graph::run_select(list<int> count, int k, float remaining_cost, float &minus_value, int step)
{
    list<int> result;
    float all_values[this->V] = {0};
    int uaf_count[this->V] = {0};
    float best_control = 0, best_impact = 0;
    for (list<int>::iterator i = count.begin(); i != count.end(); ++i)
        uaf_count[*i]++;
    float** values = new float*[this->V];
    bool** flags = new bool*[this->V];
    for (int i = 0; i < this->V; i++)
    {
        values[i] = new float[step+1];
        flags[i] = new bool[step+1];
        for (int j = 0; j < step+1; j++)
        {
            values[i][j] = 0;
            flags[i][j] = 0;
        }
        // memset(values[i], 0, sizeof(values[i]));
        // memset(flags[i], 0, sizeof(flags[i]));
    }
    for (int i = 0; i < this->V; i++)
    {
        float total_impact = this->nodes[i].p_infect * dfs(i, step, uaf_count, values, flags) * this->l_to_p;
        if (total_impact > minus_value)
            cout << i << " " << this->nodes[i].p_infect << " " << total_impact << endl;
        all_values[i] = total_impact;
        if (best_impact < total_impact)
        {
            best_control = i;
            best_impact = total_impact;
        }
    }
    bool selected[this->V] = {0};
    for (int i = 0; i < k; i++)
    {
        int to_be_added = -1;
        float best = -1;
        for (int j = 0; j < this->V; j++)
        {
            if (all_values[j] > best && selected[j] == false && remaining_cost - this->nodes[j].cost >= 0)
            {
                to_be_added = j;
                best = all_values[j];
            }
        }
        if (to_be_added == -1)
        {
            break;
        }
        result.push_back(to_be_added);
        remaining_cost -= this->nodes[to_be_added].cost;
        minus_value -= best;
        selected[to_be_added] = true;
    }
    delete[] values;
    delete[] flags;
    return result;
}

C_Ret Graph::iterative_calculate(int u, int k, bool is_infected, float p_to_l, float l_to_p)
{
    // auto a = high_resolution_clock::now();
    // static float sub = 0, sub2 = 0;
    int num_nodes = this->V;
    C_Ret result(num_nodes);
    if (k == 0)
    {
        result.g[u] = is_infected == false ? 1 : 1 - p_to_l;
        if (!is_infected)
        {
            result.h = this->nodes[u].p_infect * l_to_p;
            result.hh[u] = this->nodes[u].p_infect * l_to_p;
            result.hhh[u] = this->nodes[u].p_infect * l_to_p;
        }
        return result;
    }
    float pre_h1[num_nodes] = {0}, pre_h2[num_nodes] = {0}, pre_h3[num_nodes] = {0}, h1[num_nodes] = {0}, h2[num_nodes] = {0}, h3[num_nodes] = {0}, cur_g[num_nodes] = {0};
    pre_h2[u] = 1;
    pre_h3[u] = 1;
    result.g[u] = is_infected == false ? 1 : 1 - p_to_l;
    list<int> pre_queue, queue;
    bool queue_element[num_nodes] = {0};
    pre_queue.push_back(u);
    int cur_k = 1;
    while (true)
    {
        int cur_node = pre_queue.front();
        pre_queue.pop_front();
        if (this->adj[cur_node].size() == 0)
        {
            h1[cur_node] += pre_h1[cur_node];
            h2[cur_node] += pre_h2[cur_node];
            h3[cur_node] += pre_h3[cur_node];
            if (cur_k != k && queue_element[cur_node] == false)
            {
                queue_element[cur_node] = true;
                queue.push_back(cur_node);
            }
            if (cur_k == k && is_infected == false)
            {
                if (!this->nodes[cur_node].is_controlled)
                {
                    float tmp_p = 1 - (1 - this->nodes[cur_node].p_infect) / (1 - pre_h1[cur_node] * p_to_l);
                    result.h += pre_h1[cur_node] + pre_h2[cur_node] * (tmp_p * l_to_p);
                    result.hh[cur_node] += pre_h2[cur_node] * (tmp_p * l_to_p);
                    result.hhh[cur_node] += (pre_h1[cur_node] + pre_h2[cur_node]) * (tmp_p * l_to_p);
                }
                else{
                    result.h += pre_h1[cur_node];
                }
            }
            if (cur_k == k)
            {
                if (this->nodes[cur_node].is_controlled == true)
                    cur_g[cur_node] = 0;
                else if (is_infected == false)
                {
                    cur_g[cur_node] += pre_h1[cur_node] * p_to_l;
                }
                else if (is_infected == true)
                {
                    cur_g[cur_node] += pre_h3[cur_node] * p_to_l;
                }
            }
        }
        else{
            for (list<C_Edge>::iterator i = this->adj[cur_node].begin(); i != this->adj[cur_node].end(); ++i)
            {
                // auto cc = high_resolution_clock::now();
                int neighbor_node = i->v;
                float neighbor_weight = i->w;
                float p = this->nodes[cur_node].p_infect;
                float a = neighbor_weight * (pre_h1[cur_node] + pre_h2[cur_node] * p * l_to_p);
                float b = neighbor_weight * (pre_h2[cur_node] * (1 - p * l_to_p));
                float c = neighbor_weight *  pre_h3[cur_node];
                h1[neighbor_node] += a;
                h2[neighbor_node] += b;
                h3[neighbor_node] += c;
                if (cur_k != k && queue_element[neighbor_node] == false)
                {
                    queue.push_back(neighbor_node);
                    queue_element[neighbor_node] = true;
                }
                if (cur_k == k && is_infected == false)
                {
                    if (!this->nodes[neighbor_node].is_controlled)
                    {
                        float tmp_p = 1 - (1 - this->nodes[neighbor_node].p_infect) / (1 - a * p_to_l);
                        result.h += a + b * (tmp_p * l_to_p);
                        result.hh[neighbor_node] += b * (tmp_p * l_to_p);
                        result.hhh[neighbor_node] += (a + b) * (tmp_p * l_to_p);
                    }
                    else{
                        result.h += a;
                        result.hh[neighbor_node] = 0;
                        result.hhh[neighbor_node] = 0;
                    }
                }
                if (this->nodes[neighbor_node].is_controlled == true)
                    cur_g[neighbor_node] = 0;
                else if (is_infected == false)
                {
                    cur_g[neighbor_node] += a * p_to_l;
                }
                else if (is_infected == true)
                {
                    cur_g[neighbor_node] += c * p_to_l;
                }
                // auto dd = high_resolution_clock::now();
                // sub += (float)duration_cast<microseconds>(dd - cc).count();
            }
        }
        if (pre_queue.size() == 0)
        {
            if (cur_k == k)
                break;
            cur_k++;
            for (int i = 0; i < this->V; i++)
                result.g[i] *= 1 - cur_g[i];
            pre_queue = queue;
            queue.clear();
            memcpy(pre_h1, h1, sizeof(h1));
            memcpy(pre_h2, h2, sizeof(h2));
            memcpy(pre_h3, h3, sizeof(h3));
            memset(queue_element, 0, sizeof(queue_element));
            memset(h1, 0, sizeof(h1));
            memset(h2, 0, sizeof(h2));
            memset(h3, 0, sizeof(h3));
            memset(cur_g, 0, sizeof(cur_g));
        }  
    }
    // auto b = high_resolution_clock::now();
    // sub2 += (float)duration_cast<microseconds>(b - a).count();
    // cout << sub << " " << sub2 << " " << sub / sub2 << endl;
    return result;
}
// Driver program to test above 
void g_main(int num_affected, int num_maxs, int step, float p_to_l, float l_to_p, float cost_control, string graph_name, string uaf_distribution_name, string cost_path, bool lower_bound, int num_per_iter, string prefix, bool ls)
{
    map<string, bool> nodes;
    string name = "./" + graph_name + ".txt";
    ifstream f;
    f.open(name.c_str());
    string line;
    while (getline(f, line))
    {
        int pos = line.find('\t');
        string u = line.substr(0, pos);
        int pos2 = line.substr(pos+1, line.length()).find('\t');
        string v = line.substr(pos+1, line.length()).substr(0, pos2);
        float w = stof(line.substr(pos+1, line.length()).substr(pos2+1, line.substr(pos+1, line.length()).length()));
        nodes[u] = true;
        nodes[v] = true;
    }
    f.close();
    cout << "Total number of POIs: " << nodes.size() <<endl;
    ff << "Total number of POIs: " << nodes.size() <<endl;

    list<int> uaf;
    fstream f_uninfect(uaf_distribution_name);
    while (getline(f_uninfect, line))
    {
        int u = stoi(line);
        uaf.push_back(u);
    }
    int num_unaffected = uaf.size();
    cout << "uaf: " << num_unaffected << endl;  //xiao
    Graph g(nodes.size(), num_affected, num_unaffected, p_to_l, l_to_p);

    int uaf_count[g.V] = {0};
    for (list<int>::iterator i = uaf.begin(); i != uaf.end(); ++i)
    {
        uaf_count[*i]++;
    }

    f.open(name.c_str());
    while (getline(f, line))
    {
        int pos = line.find('\t');
        int u = stoi(line.substr(0, pos));
        int pos2 = line.substr(pos+1, line.length()).find('\t');
        int v = stoi(line.substr(pos+1, line.length()).substr(0, pos2));
        float w = stof(line.substr(pos+1, line.length()).substr(pos2+1, line.substr(pos+1, line.length()).length()));
        g.addEdge(u, v, w); 
    }
    f.close();

    g.nodes = new C_Node[g.V];

    f.open(cost_path);
    while (getline(f, line))
    {
        int pos = line.find('\t');
        int u = stoi(line.substr(0, pos));
        float c = stof(line.substr(pos+1, line.find('\n')));
        g.nodes[u].cost = c;
    }
    f.close();

    float tmp_cost = 0;
    map<int, bool> current_controlled;

    auto head_timer = high_resolution_clock::now();

    while (true)
    {
        auto prev_timer = high_resolution_clock::now();
        static int itr = 0;
        cout << "~~~~~New Iteration " << itr++ << "~~~~" << endl;
        ff << "~~~~~New Iteration " << itr << "~~~~" << endl;
        float max_people = 0;

        for (map<int, bool>::iterator i = current_controlled.begin(); i != current_controlled.end(); ++i)
        {
            g.nodes[i->first].is_controlled = true;
        }
        cout << "Total number of controlled POIs: " << current_controlled.size() << endl;
        ff << "Total number of controlled POIs: " << current_controlled.size() << endl;
        list<int> best_control;
        vector<float> best_hh;
        vector<float> best_hhh;
        for (int i = 0; i < g.V; i++)
        {
            best_hh.push_back(0);
            best_hhh.push_back(0);
        }

        for (int INNER_LOOP = 0; INNER_LOOP < num_maxs; INNER_LOOP++)
        {
            cout << "max: " << INNER_LOOP << endl;
            ff << "max: " << INNER_LOOP << endl;
            float exp_nodes = 0;
            float exp_people = -10000;
            list<int> af;
            ifstream f;
            // string prefix = "../Data/";
            // string prefix = "";
            if (INNER_LOOP == 0)
                f.open(prefix + "Degree_Ranked.txt");
            if (INNER_LOOP == 1)
                f.open(prefix + "Cost_Ranked.txt");
            if (INNER_LOOP == 2)
                f.open(prefix + "Random_Ranked.txt");
            string line;
            while (getline(f, line))
            {
                if (af.size() >= num_affected)
                    break;
                af.push_back(stoi(line.substr(0, line.find('\n'))));
            }
            f.close();
            if (INNER_LOOP == 3)
            {
                for (int i = 0; i < g.V; i++)
                    af.push_back(i);
            }
            valarray<float> af_g(1, g.V);
            if (ls)
            {
                for (int i = 0; i < g.V; i++)
                {
                    g.nodes[i].p_infect = 0;
                }
            }
            float af_count[g.V] = {0};
            for (list<int>::iterator i = af.begin(); i != af.end(); ++i)
            {
                af_count[*i]++;
                if (INNER_LOOP == 3)
                    af_count[*i] = (float)num_affected / (float)g.V;
            }
            C_Ret ret;
            #pragma omp declare reduction(*:valarray<float>:omp_out *= omp_in) initializer (omp_priv = omp_orig)
            #pragma omp declare reduction(vec_float_plus : std::vector<float> : \
                          std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<float>())) \
                initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
            #pragma omp parallel for private(ret)
            for (int i = 0; i < g.V; i++)
            {
                if (af_count[i] != 0)
                {
                    C_Ret ret = g.iterative_calculate(i, step, true, p_to_l, l_to_p);
                    af_g = af_g * pow(valarray<float>(ret.g, g.V), valarray<float>(af_count[i], g.V));
                }
            }
            int inner_count = 0;
            // vector<float> overall_hh;
            int run_times = 0;
            while (true)
            {
                static int all_count = 0;
                float cur_people = 0;
                valarray<float> tmp_g(af_g);
                vector<float> tmp_hh;
                vector<float> tmp_hhh;
                for (int i = 0; i < g.V; i++)
                {
                    tmp_hh.push_back(0);
                    tmp_hhh.push_back(0);
                }
                // pthread_t threads[g.V];
                // struct thread_data td;
                // td.g = &g;
                // td.k = step;
                // td.is_infected = false;
                // td.p_to_l = p_to_l;
                // td.l_to_p = l_to_p;
                // td.cur_people = &cur_people;
                // td.tmp_g = &tmp_g;
                // td.tmp_hh = &tmp_hh;

                C_Ret ret;
                int j;
                // omp_set_dynamic(0);
                // omp_set_num_threads(8);
                #pragma omp declare reduction(*:valarray<float>:omp_out *= omp_in) initializer (omp_priv = omp_orig)
                #pragma omp declare reduction(vec_float_plus : std::vector<float> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<float>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
                #pragma omp parallel for private(ret, j) reduction(+:cur_people)
                    for (int i = 0; i < g.V; i++)
                    {
                        // cout << i << endl;
                        if (uaf_count[i] != 0)
                        {
                            // td.u = i;
                            // td.uaf_count = uaf_count[i];        
                            // pthread_create(&threads[i], NULL, Func, (void *)&td);
                            // Func((void *) &td);
                            C_Ret ret = g.iterative_calculate(i, step, false, p_to_l, l_to_p);
                            cur_people += uaf_count[i] * ret.h;
                            tmp_g = tmp_g * pow(valarray<float>(ret.g, g.V), valarray<float>(uaf_count[i], g.V));
                            for (j = 0; j < g.V; j++)
                            {
                                tmp_hh[j] += uaf_count[i] * ret.hh[j];
                                tmp_hhh[j] += uaf_count[i] * ret.hhh[j];
                            }
                            // delete[] ret.g;
                            // delete[] ret.hh;
                        }
                    }
                inner_count++;
                all_count++;
                // cout << cur_people << "\t" << exp_people << "\t" << inner_count << "\t" << all_count << endl;
                // ff << cur_people << "\t" << exp_people << "\t" << inner_count << "\t" << all_count << endl;
                if ((abs(exp_people) < 1e-6 && abs(cur_people) < 1e-6) || abs(cur_people - exp_people)/cur_people < 1e-4 || (ls && run_times++ > 0))
                {
                    cout << cur_people << "\t" << exp_people << "\t" << inner_count << "\t" << all_count << endl;
                    ff << cur_people << "\t" << exp_people << "\t" << inner_count << "\t" << all_count << endl;
                    exp_people = cur_people;
                    for (int i = 0; i < g.V; i++)
                    {
                        if (best_hh[i] < cur_people - tmp_hh[i])
                            best_hh[i] = cur_people - tmp_hh[i];
                        if (best_hhh[i] < cur_people - tmp_hhh[i])
                            best_hhh[i] = cur_people - tmp_hhh[i];
                    }
                    break;
                }
                exp_people = cur_people;
                for (int i = 0; i < g.V; i++)
                {
                    if (g.nodes[i].is_controlled == false)
                    {
                        g.nodes[i].p_infect = 1 - tmp_g[i];
                        // if (INNER_LOOP == 3)
                        // {
                        //     g.nodes[i].p_infect += (float)g.num_affected / (float)g.V;
                        // }
                        if (g.nodes[i].p_infect < 1e-16)
                            g.nodes[i].p_infect = 0;
                    }
                    else
                    {
                        g.nodes[i].p_infect = 0;
                    }
                }
            }

            // float minus_value = exp_people;
            // list<int> select_set;

            // if (lower_bound)
            // {
            //     float inner_cost = 0;
            //     for (auto i: sort_indexes(overall_hh))
            //     {
            //         if (tmp_cost + g.nodes[i].cost + inner_cost > cost_control)
            //             continue;
            //         select_set.push_back(i);
            //         inner_cost += g.nodes[i].cost;
            //         minus_value -= overall_hh[i];
            //         // cout << i << " " << overall_hh[i] << endl;
            //         if (select_set.size() >= num_per_iter)
            //             break;
            //     }
            // }else{
            //     select_set = g.run_select(uaf, num_per_iter, cost_control - tmp_cost, minus_value, step);
            // }
            // if (minus_value == 0 && exp_people > 0)
            // {
            //     best_control = select_set;
            // }
            // cout << "Temp max: " << INNER_LOOP << " " << max_people << " " << exp_people << " " << minus_value << " " << exp_people - minus_value << endl;
            // if (max_people < exp_people)
            // {
            //     max_people = exp_people;
            //     best_control = select_set;
            //     // cout << "Temp max: " << INNER_LOOP << " " << exp_people << endl;
            //     ff << "Temp max: " << INNER_LOOP << " " << exp_people << endl;
            // }
        }
        list<int> select_set;
        int ccc = 0;
        for (auto i: sort_indexes(best_hh))
        {
            cout << best_hh[i] << " " << i << endl;
            if (ccc++ > 5)
                break;
        }
        if (lower_bound)
        {
            float inner_cost = 0;
            for (auto i: sort_indexes(best_hh))
            {
                if (best_hh[i] < 1e-6 || g.nodes[i].is_controlled)
                    continue;
                if (tmp_cost + g.nodes[i].cost + inner_cost > cost_control)
                    continue;
                select_set.push_back(i);
                inner_cost += g.nodes[i].cost;
                // cout << i << " " << overall_hh[i] << endl;
                if (select_set.size() >= num_per_iter)
                    break;
            }
        }else{
            float inner_cost = 0;
            for (auto i: sort_indexes(best_hhh))
            {
                if (best_hhh[i] < 1e-6 || g.nodes[i].is_controlled)
                    continue;
                if (tmp_cost + g.nodes[i].cost + inner_cost > cost_control)
                    continue;
                select_set.push_back(i);
                inner_cost += g.nodes[i].cost;
                // cout << i << " " << overall_hh[i] << endl;
                if (select_set.size() >= num_per_iter)
                    break;
            }
        }
        best_control = select_set;
        if (best_control.size() == 0)
            break;
        for (list<int>::iterator i = best_control.begin(); i != best_control.end(); ++i)
        {
            cout << "Adding: " << *i << endl;
            ff << "Adding: " << *i << endl;
            current_controlled[*i] = true;
            tmp_cost += g.nodes[*i].cost;

        }
        auto now_timer = high_resolution_clock::now();
        cout << "Total time for this iteration: " << (float)duration_cast<microseconds>(now_timer - prev_timer).count()/1e6 << " " << (float)duration_cast<microseconds>(now_timer - head_timer).count()/1e6 << endl;
        ff   << "Total time for this iteration: " << (float)duration_cast<microseconds>(now_timer - prev_timer).count()/1e6 << " " << (float)duration_cast<microseconds>(now_timer - head_timer).count()/1e6 << endl;
    }
    for (map<int, bool>::iterator i = current_controlled.begin(); i != current_controlled.end(); ++i)
    {
        ff << i->first << endl;
        cout << i->first << endl;
    }

} 

int main(int argc, char* argv[])
{

    srand((unsigned) time(0));

    // float cost_control = 500000;
    // int num_affected = 50;
    // int num_per_iter = 10000;
    // // int num_unaffected = 20000;
    // int num_unaffected = -1;
    // int step = 3;
    // float p_to_l = 0.1;
    // float l_to_p = 0.1;
    // ff.open("log_true_all.txt");
    // string prefix = "../Data/";
    // bool lower_bound = false;

    float cost_control = atof(argv[1]);
    int num_affected = atoi(argv[2]);
    int num_per_iter = atoi(argv[3]);
    int step = atoi(argv[4]);
    float p_to_l = atof(argv[5]);
    float l_to_p = atof(argv[6]);
    ff.open(string(argv[7]));
    string prefix = string(argv[8]);
    bool lower_bound = atoi(argv[9]);
    bool ls = atoi(argv[10]);

    string uaf_mult = string(argv[11]);
    int num_maxs = 4;     // initializations of infected/uninfected people changed
    string graph_name = prefix + "Edges";
    // string graph_name = "graph_SCC_weighted";
    string uaf_distribution_name = prefix + uaf_mult + "Spectator_Distribution.txt";
    string cost_path = prefix + "Cost.txt";

    g_main(num_affected, num_maxs, step, p_to_l, l_to_p, cost_control, graph_name, uaf_distribution_name, cost_path, lower_bound, num_per_iter, prefix, ls);

    ff.close();
    return 0;
}
