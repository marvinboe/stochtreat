/*
 *  dependency.cpp
 *  stochmut
 *
 *  Created by Tom Lenaerts on 1/31/11.
 *  Copyright 2011 Universite Libre de Bruxelles. All rights reserved.
 *
 */

#include "dependency.h"

DependencyNode::DependencyNode(const DependencyNode& other){
    _affecting.clear();
    _idx = other.reaction();
    std::vector<DependencyNode*>::const_iterator start = other.begin();
    std::vector<DependencyNode*>::const_iterator stop = other.end();
    while(start !=stop){
        affects(*start);
        ++start;
    }
}

DependencyNode& DependencyNode::operator=(const DependencyNode& other){
    _affecting.clear();
    _idx =other.reaction();
    std::vector<DependencyNode*>::const_iterator start = other.begin();
    std::vector<DependencyNode*>::const_iterator stop = other.end();
    while(start !=stop){
        affects(*start);
        ++start;
    }
    return *this;
}

std::ostream& DependencyNode::display(std::ostream& os){
    os << "<" << _idx << "|"<< _affecting.size()<<"|";
    for(unsigned i=0; i< _affecting.size(); ++i){
        unsigned tmp =  _affecting[i]->reaction();
        os << tmp;
        if(i < (_affecting.size() -1))
            os << ",";
    }
    os << ">";
    return os;
}



//need to add a mutation reaction at a later stage; mutations occur between cell types. 
DependencyGraph::DependencyGraph(Model& pool, Data& data, AllReactions& all, RanGen& ran){
    //for every type of cell in every compartment there is one self-renewal and one differentiation reaction
    //note that compartment 0 is a special case, since it's size needs to remain constant
    unsigned pos;
    double sum=0.0;
    unsigned numtypes = 4;
    for(unsigned k=1; k < pool.numStoch(); ++k) { //every compartment

        //
        Treatment *treat= new Treatment(k,pool.getTreatRate());// type(k) -> type(k+1)+type(k+1)
        treat->setPropensity((treat->sufficientReactants(pool)?treat->reactantFactor(pool):0.0)); 
        sum += treat->propensity();
        pos=all.add(treat);
        DependencyNode *treatnode=new DependencyNode(pos);
        treatnode->affects(treatnode);

        for(unsigned cell_type_id = 0; cell_type_id < numtypes; ++cell_type_id) {
            // cout << "(1-eps) = " << (1.0-data.eps(cell_type_id)) << "\t rate= " << pool.getRate(k,cell_type_id) << endl; 
            SelfRenewal *self= new SelfRenewal(k,cell_type_id,(1.0 - data.eps(cell_type_id))*pool.getRate(k,cell_type_id));// cell_type_id(k) -> cell_type_id(k)+cell_type_id(k)
            if (k<=data.n_neutral_compartments()){
                self->setRate((1.0 - data.eps(0))*pool.getRate(k,0));
            }
            self->setPropensity((self->sufficientReactants(pool)?self->reactantFactor(pool):0.0)); 
            sum += self->propensity();
            pos=all.add(self);
            DependencyNode *selfnode=new DependencyNode(pos);

            Differentation *diff= new Differentation(k,cell_type_id,data.eps(cell_type_id)*pool.getRate(k,cell_type_id));// cell_type_id(k) -> cell_type_id(k+1)+cell_type_id(k+1)
            if (k<=data.n_neutral_compartments()){
                diff->setRate(data.eps(0)*pool.getRate(k,0));
            }
            diff->setPropensity((diff->sufficientReactants(pool)?diff->reactantFactor(pool):0.0)); 
            sum += diff->propensity();
            pos=all.add(diff);
            DependencyNode *diffnode=new DependencyNode(pos);


            //add dependencies between both nodes and themselves
            selfnode->affects(selfnode);
            selfnode->affects(diffnode);
            diffnode->affects(diffnode);
            diffnode->affects(selfnode);
            if (cell_type_id==1 || cell_type_id==3){
                treatnode->affects(selfnode);
                treatnode->affects(diffnode);
                if (cell_type_id==1) { // cell_type_id is cancer cell
                 selfnode->affects(treatnode);
                 diffnode->affects(treatnode);
                }
            }


            //both nodes are also affected by certain reactions from the previous compartment.			
            //need to add them. Carefull with the reactions in compartment 0.
            if( (k-1) > 0){ //starting from k=2
                _diffnodes[((k-2)*4 + cell_type_id)]->affects(selfnode);
                _diffnodes[((k-2)*4 + cell_type_id)]->affects(diffnode);
                if (cell_type_id==1){
                    _diffnodes[((k-2)*4 + cell_type_id)]->affects(treatnode);
                }
            }
            //add both nodes to their corresponding vectors.
            pos=add(SELF_RENEWAL,selfnode);
            self->setDG(SELF_RENEWAL,pos);
            pos=add(DIFFERENTATION,diffnode);
            diff->setDG(DIFFERENTATION, pos);	
        }

        pos=add(TREATMENT,treatnode);//add treatment node to std::vector
        treat->setDG(TREATMENT,pos);
    }
    //add reactions for the stem cell compartment:
    for(int first_tp = 0; first_tp < 3; ++first_tp) { // excluding the bound type, i.e.tp = 3
        //differentiation with renewal (->Moran step)
        StemCellRenewal *scr=new StemCellRenewal(&ran,first_tp,
                data.beta()*pool.getRate(0,first_tp));
        scr->setPropensity(scr->sufficientReactants(pool)?scr->reactantFactor(pool):0.0);
        sum += scr->propensity();
        pos=all.add(scr);
        DependencyNode *scrnode=new DependencyNode(pos);
        pos=add(MORAN,scrnode);
        scr->setDG(MORAN, pos);

        //		cout << "# stem cell renewal : " << *scr << "\t propensity " << <<endl;
        scrnode->affects(scrnode);

        //differentiation without renewal
        StemCellDifferentiation *scd=new StemCellDifferentiation(first_tp,(1-data.beta())*pool.getRate(0,first_tp));
        scd->setPropensity(scd->sufficientReactants(pool)?scd->reactantFactor(pool):0.0);
        sum += scd->propensity();
        pos=all.add(scd);
        DependencyNode *scdnode=new DependencyNode(pos);
        pos=add(MORAN,scdnode);
        scd->setDG(MORAN, pos);

        //		cout << "# stem cell renewal : " << *scr << "\t propensity " << <<endl;
        scdnode->affects(scdnode);

        if (_selfnodes.size() > 0 && _diffnodes.size() > 0){
            //Stem Cell reaction moves random cell to next comp. -> adding links to all other cell types
            DependencyNode *nexttreat=get(TREATMENT,0);
            scrnode->affects(nexttreat);
            scdnode->affects(nexttreat);
            for(int second_type = 0; second_type < 3; ++second_type) {
                DependencyNode *nextself, *nextdiff;
                nextself = get(SELF_RENEWAL,second_type); 
                nextdiff = get(DIFFERENTATION,second_type);
                scrnode->affects(nextself);
                //				cout << *(all[scrnode->reaction()]) << "is linked to " << *(all[nextself->reaction()]) << endl;
                scrnode->affects(nextdiff);
                //				cout << *(all[scrnode->reaction()]) << "is linked to " << *(all[nextdiff->reaction()]) << endl;
                scdnode->affects(nextself);
                scdnode->affects(nextdiff);
            }
        }

    }

    //add dependencies between all reactions in first (Moran) compartment
    for(unsigned first_tp = 0; first_tp < 3; ++first_tp) { 
        DependencyNode *scrnode =  get(MORAN,first_tp*2);
        for(unsigned second_tp = 0; second_tp < 3; ++second_tp) {
            if(first_tp != second_tp){
                DependencyNode *oscrnode =  get(MORAN,second_tp*2);
                scrnode->affects(oscrnode);
                //				cout << *(all[scrnode->reaction()]) << "is linked to " << *(all[oscrnode->reaction()]) << endl;

            }
            DependencyNode *scdnode =  get(MORAN,second_tp*2+1);
            scrnode->affects(scdnode);
        }
    }
    all.setPropSum(sum);
    // display(std::cout);
    // cout << "#finished creation DependencyGraph, number of reactions = "<< all.size() << endl;
}


std::ostream& DependencyGraph::display(std::ostream& os){
    os << "DependencyGraph [" << std::endl << "Moran reactions";
    std::vector<DependencyNode*>::iterator startmoran = _morannodes.begin();
    std::vector<DependencyNode*>::iterator stopmoran = _morannodes.end();
    while(startmoran !=stopmoran){
        os << std::endl << **startmoran;
        ++startmoran;
    }
    os << std::endl<< "Selfrenewal";
    std::vector<DependencyNode*>::iterator startself = _selfnodes.begin();
    std::vector<DependencyNode*>::iterator stopself = _selfnodes.end();
    while(startself !=stopself){
        os << std::endl << **startself;
        ++startself;
    }
    os << std::endl << "Differentiation";
    std::vector<DependencyNode*>::iterator startdiff = _diffnodes.begin();
    std::vector<DependencyNode*>::iterator stopdiff = _diffnodes.end();
    while(startdiff !=stopdiff){
        os << std::endl << **startdiff;
        ++startdiff;
    }

    os << std::endl << "Treatment";
    std::vector<DependencyNode*>::iterator treat_it = _treatnodes.begin();
    std::vector<DependencyNode*>::iterator stoptreat = _treatnodes.end();
    while(treat_it !=stoptreat){
        os << std::endl << **treat_it;
        ++treat_it;
    }
    os << "]"<<std::endl;
    return os;
}




