/*
 *  dependency.h
 *  stochmut
 *
 *  Created by Tom Lenaerts on 1/31/11.
 *  Copyright 2011 Universite Libre de Bruxelles. All rights reserved.
 *
 */

#ifndef __DEPENDENCY_H
#define __DEPENDENCY_H

#include <string>
#include <istream>
#include <iostream>
#include <vector>
#include "model.h"
#include "data.h"
#include "reactions.h"


class DependencyNode {
    public:
        DependencyNode(){_affecting.clear();}
        DependencyNode(unsigned idx):_idx(idx){_affecting.clear();}
        DependencyNode(const DependencyNode& other);
        ~DependencyNode(){
            _affecting.clear();
        }
        DependencyNode& operator=(const DependencyNode& other);

        /** returns index of reaction in list of all reactions. */
        unsigned reaction() const {return _idx;}

        /** adds pointer to other DependencyNode to _affecting std::vector.*/
        void affects(DependencyNode* other) {
            _affecting.push_back(other);
        }

        std::vector<DependencyNode*>::iterator begin() {return _affecting.begin();}
        std::vector<DependencyNode*>::iterator end() {return _affecting.end();}
        std::vector<DependencyNode*>::const_iterator begin() const {return _affecting.begin();}
        std::vector<DependencyNode*>::const_iterator end() const {return _affecting.end();}
        friend std::ostream & operator<<(std::ostream &o, DependencyNode& dn){return dn.display(o);}

    protected:
        std::ostream& display(std::ostream& os);
        unsigned _idx; /** < index in list of all reactions. */
        std::vector<DependencyNode*> _affecting;
};

enum reaction_type {MORAN=0,SELF_RENEWAL,DIFFERENTATION,TREATMENT};

class DependencyGraph {
    public:
        DependencyGraph(){_morannodes.clear(); _selfnodes.clear(); _diffnodes.clear();}
        DependencyGraph(Model&, Data&, AllReactions&, RanGen& ran);
        //	DependencyGraph(Model&, Data&, AllReactions&);
        ~DependencyGraph(){
            while(_morannodes.size() > 0){
                DependencyNode* dn = _morannodes[_morannodes.size()-1];
                _morannodes.pop_back();
                delete dn;
            }
            while(_selfnodes.size() > 0){
                DependencyNode* dn = _selfnodes[_selfnodes.size()-1];
                _selfnodes.pop_back();
                delete dn;
            }
            while(_diffnodes.size() > 0){
                DependencyNode* dn = _diffnodes[_diffnodes.size()-1];
                _diffnodes.pop_back();
                delete dn;
            }
            while(_treatnodes.size() > 0){
                DependencyNode* dn = _treatnodes[_treatnodes.size()-1];
                _treatnodes.pop_back();
                delete dn;
            }
        }

        unsigned numMoranNodes() const {return (unsigned) _morannodes.size();}
        unsigned numSelfNodes() const {return (unsigned) _selfnodes.size();}
        unsigned numDiffNodes() const {return (unsigned) _diffnodes.size();}

        DependencyNode* getMoranNode(unsigned pos) {return _morannodes[pos];}
        DependencyNode* getSelfNode(unsigned pos) {return _selfnodes[pos];}
        DependencyNode* getDiffNode(unsigned pos) {return _diffnodes[pos];}
        DependencyNode* getTreatNode(unsigned pos) {return _treatnodes[pos];}

        /** Returns the dependency node object corresponding to the 
         * type of the reaction and the number in the reaction-specific
         * container. */
        DependencyNode* get(unsigned which, unsigned pos){
            switch(which){
                case MORAN:
                    return getMoranNode(pos); 
                case SELF_RENEWAL: 
                    return getSelfNode(pos);
                case DIFFERENTATION:
                    return getDiffNode(pos);
                case TREATMENT:
                    return getTreatNode(pos);
            }
            return NULL;
        }

        unsigned addMoranNode(DependencyNode* dn) {
            _morannodes.push_back(dn);
            return (unsigned) _morannodes.size()-1;
        }
        unsigned addSelfNode(DependencyNode* dn) {
            _selfnodes.push_back(dn);
            return (unsigned) _selfnodes.size()-1;
        }
        unsigned addDiffNode(DependencyNode* dn) {
            _diffnodes.push_back(dn);
            return  (unsigned) _diffnodes.size()-1;
        }

        unsigned addTreatNode(DependencyNode* dn) {
            _treatnodes.push_back(dn);
            return  (unsigned) _treatnodes.size()-1;
        }

        unsigned add(int which, DependencyNode* dn){
            switch(which){
                case MORAN:
                    return addMoranNode(dn); 
                case SELF_RENEWAL: 
                    return addSelfNode(dn);
                case DIFFERENTATION:
                    return addDiffNode(dn);
                case TREATMENT:
                    return addTreatNode(dn);
            }
            return 999;
        }

        friend std::ostream & operator<<(std::ostream &o, DependencyGraph& dg){return dg.display(o);}

    private:
        std::ostream& display(std::ostream& os);
        std::vector<DependencyNode*> _morannodes;
        std::vector<DependencyNode*> _selfnodes;
        std::vector<DependencyNode*> _diffnodes;
        std::vector<DependencyNode*> _treatnodes;
};


#endif
