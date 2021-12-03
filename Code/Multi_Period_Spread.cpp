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
// #include <omp.h>
using namespace std; 
#include <chrono> 
#include <cmath>
using namespace std::chrono; 
  
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
    C_Ret(int V)
    {
        h = 0;
        g = new float[V];
        for (int i = 0; i < V; i++)
            this->g[i] = 1;
    }
    C_Ret()
    {
        g = new float[2];
    }
    ~C_Ret()
    {
        delete[] g;
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
    C_Ret iterative_calculate(int u, int k, bool is_infected, float p_to_l, float l_to_p, float initial_s);
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

C_Ret Graph::iterative_calculate(int u, int k, bool is_infected, float p_to_l, float l_to_p, float initial_s)
{
    int num_nodes = this->V;
    C_Ret result(num_nodes);
    if (k == 0)
    {
        result.g[u] = is_infected == false ? 1 : 1 - p_to_l;
        if (!is_infected)
        {
            result.h = this->nodes[u].p_infect * l_to_p;
        }
        return result;
    }
    // float pre_h1[num_nodes] = {0}, pre_h2[num_nodes] = {0}, pre_h3[num_nodes] = {0}, h1[num_nodes] = {0}, h2[num_nodes] = {0}, h3[num_nodes] = {0}, cur_g[num_nodes] = {0};
    vector<float> pre_h1(num_nodes, 0), pre_h2(num_nodes, 0), pre_h3(num_nodes, 0), h1(num_nodes, 0), h2(num_nodes, 0), h3(num_nodes, 0), cur_g(num_nodes, 0);
    pre_h1[u] = initial_s;
    pre_h2[u] = 1 - initial_s;
    pre_h3[u] = 1;
    result.g[u] = is_infected == false ? 1 : 1 - p_to_l;
    list<int> pre_queue, queue;
    // bool queue_element[num_nodes] = {0};
    vector<bool> queue_element(num_nodes, 0);
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
                    }
                    else{
                        result.h += a;
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
            // memcpy(pre_h1, h1, sizeof(h1));
            // memcpy(pre_h2, h2, sizeof(h2));
            // memcpy(pre_h3, h3, sizeof(h3));
            pre_h1 = h1;
            pre_h2 = h2;
            pre_h3 = h3;
            fill(queue_element.begin(), queue_element.end(), 0);
            fill(h1.begin(), h1.end(), 0);
            fill(h2.begin(), h2.end(), 0);
            fill(h3.begin(), h3.end(), 0);
            fill(cur_g.begin(), cur_g.end(), 0);
            // memset(queue_element, 0, sizeof(queue_element));
            // memset(h1, 0, sizeof(h1));
            // memset(h2, 0, sizeof(h2));
            // memset(h3, 0, sizeof(h3));
            // memset(cur_g, 0, sizeof(cur_g));
        }  
    }
    return result;
}
// Driver program to test above 
void g_main(int num_affected, int num_maxs, int step, float p_to_l, float l_to_p, float cost_control, string c_path, string graph_name, string uaf_distribution_name, string cost_path, string prefix)
{ 
    ff << c_path << endl;
    cout << c_path << endl;
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
    cout << "uaf: " << num_unaffected << endl; //xiao
    Graph g(nodes.size(), num_affected, num_unaffected, p_to_l, l_to_p);

    // int uaf_count[g.V] = {0};
    vector<int> uaf_count(g.V, 0);
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


    float max_people = 0;
    float max_location = 0;

    // float location_initial_s[g.V] = {0};
    list<int> af;
    f.open(prefix + "Degree_Ranked.txt");
    while (getline(f, line))
    {
        if (af.size() >= num_affected)
            break;
        af.push_back(stoi(line.substr(0, line.find('\n'))));
    }
    f.close();
    
    // float af_count[g.V] = {0};
    vector<float> af_count(g.V, 0);
    for (list<int>::iterator i = af.begin(); i != af.end(); ++i)
    {
        af_count[*i]++;
    }

    int incubation_period = 14;
    vector< vector<float> > global_infected(g.V, vector<float>(incubation_period, 0));
    int head = 0, tail = 0;
    for (int i = 0; i < g.V; i++)
        global_infected[i][tail] = af_count[i];
    tail++;
    cout << head << " " << tail << endl;

    // Start multi-period spread process
    vector<int> current_controlled;
    for (int INNER_LOOP = 0; INNER_LOOP < 180; INNER_LOOP++)
    {
        current_controlled.clear();

        cout << "Here1" << endl;
        float tmp_cost = 0;
        ifstream f;
        string line;
        // cout << "Here" << endl;
        f.open(c_path);
        // cout << "Here1" << endl;
        while (getline(f, line))
        {
            int u = stoi(line.substr(0, line.find('\n')));
            // cout << tmp_cost << " " << g.nodes[u].cost << " "<< cost_control << endl;
            if (tmp_cost + g.nodes[u].cost <= cost_control)
            {
                // cout << "Here2" << endl;
                current_controlled.push_back(u);
                tmp_cost += g.nodes[u].cost;
            }
            
        }
        f.close();

        for (int i = 0; i < g.V; i++)
        {
            g.nodes[i].is_controlled = false;
        }

        for (auto i = current_controlled.begin(); i != current_controlled.end(); ++i)
        {
            int ii = *i;
            g.nodes[ii].is_controlled = true;
        }
        cout << "Total number of controlled POIs: " << current_controlled.size() << endl;
        ff << "Total number of controlled POIs: " << current_controlled.size() << endl;

        C_Ret ret_inf;
        valarray<float> af_g(1, g.V);
        // cout << "Here" << endl;
        // #pragma omp declare reduction(*:valarray<float>:omp_out *= omp_in) initializer (omp_priv = omp_orig)
        // #pragma omp declare reduction(vec_float_plus : std::vector<float> : \
        //               std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<float>())) \
        //     initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
        // #pragma omp parallel for private(ret_inf)
        for (int i = 0; i < g.V; i++)
        {
            if (af_count[i] != 0)
            {
                C_Ret ret_inf = g.iterative_calculate(i, step, true, p_to_l, l_to_p, 1);
                af_g = af_g * pow(valarray<float>(ret_inf.g, g.V), valarray<float>(af_count[i], g.V));
            }
        }
        // cout << "Here" << endl;
        float exp_nodes = 0;
        float exp_people = -10000;
        float exp_location = 0;
        int inner_count = 0;
        // for (int i = 0; i < 100; i++)
        //     cout << location_initial_s[i] << " ";
        // cout << endl;
        while (true)
        {
            // float location_current_s[g.V] = {0};
            vector<float> location_current_s(g.V, 0);
            static int all_count = 0;
            
            float cur_people = 0;
            float cur_location = 0;
            valarray<float> tmp_g(af_g);
            // if (INNER_LOOP == 0)
            //     tmp_g = new valarray<float>(af_g);
            // else
            //     tmp_g = new valarray<float>(1, g.V);
            C_Ret ret;
            // #pragma omp declare reduction(*:valarray<float>:omp_out *= omp_in) initializer (omp_priv = omp_orig)
            // #pragma omp declare reduction(vec_float_plus : std::vector<float> : \
            //               std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<float>())) \
            //     initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
            // #pragma omp parallel for private(ret) reduction(+:cur_people)
            for (int i = 0; i < g.V; i++)
            {
                if (uaf_count[i] != 0)
                {
                    C_Ret ret = g.iterative_calculate(i, step, false, p_to_l, l_to_p, 0);
                    location_current_s[i] = ret.h;
                    cur_people += uaf_count[i] * (ret.h - 0);
                    tmp_g = tmp_g * pow(valarray<float>(ret.g, g.V), valarray<float>(uaf_count[i], g.V));
                }
            }
            cout << cur_people << "\t" << exp_people << "\t" << inner_count << "\t" << all_count++ << endl;
            ff << cur_people << "\t" << exp_people << "\t" << inner_count++ << "\t" << all_count << endl;
            if ((abs(exp_people) < 1e-6 && abs(cur_people) < 1e-6) || abs(cur_people - exp_people)/cur_people < 1e-4)
            {
                // cout << "Day " << INNER_LOOP << ": " << cur_people << endl;
                exp_people = 0;
                for (int i = 0; i < g.V; i++)
                {
                    if (g.nodes[i].is_controlled == false)
                        cur_location += 1 - tmp_g[i];
                }
                exp_location = cur_location;

                for (int i = 0; i < g.V; i++)
                {
                    float location_new_infected = uaf_count[i] * location_current_s[i];
                    uaf_count[i] -= location_new_infected;
                    af_count[i] += location_new_infected;
                    if (tail == head)
                    {
                        af_count[i] -= global_infected[i][head];
                        exp_people += global_infected[i][head];
                    }
                    global_infected[i][tail] = location_new_infected;
                }
                if (head == tail)
                    head = (head+1) % incubation_period;
                tail = (tail+1) % incubation_period;
                cout << "here" << endl;
                float w = 0.01, h = -4;
                cost_control = 1.0/(1 + exp(-1 * (w * exp_people + h)));
                cout << "Day " << INNER_LOOP << " newly infected: " << exp_people << " " << head << " " << tail << ", control budget: " << cost_control << endl;
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
        cout << "here" << endl;
        // if (max_people < exp_people)
        // {
        //     max_people = exp_people;
        //     max_location = exp_location;
        //     cout << "Temp max: " << INNER_LOOP << " " << exp_people << " " << exp_location << endl;
        //     ff << "Temp max: " << INNER_LOOP << " " << exp_people << " " << exp_location << endl;
        // }
    }
    ff << "Final result: " << max_people << " " << max_location << endl;
    cout << "Final result: " << max_people << " " << max_location << endl;
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
    string control_strategy_file_name = prefix + "Cost_Ranked.txt";

    string uaf_mult = string(argv[11]);
    int num_maxs = 4;     // initializations of infected/uninfected people changed
    string graph_name = prefix + "Edges";
    // string graph_name = "graph_SCC_weighted";
    string uaf_distribution_name = prefix + uaf_mult + "Spectator_Distribution.txt";
    string cost_path = prefix + "Cost.txt";

    g_main(num_affected, num_maxs, step, p_to_l, l_to_p, cost_control, control_strategy_file_name, graph_name, uaf_distribution_name, cost_path, prefix);

    ff.close();
    return 0;
}
