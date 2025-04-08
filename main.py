from typing import Union 
from fastapi import FastAPI
from minesweeper import MinesweeperSolver
from pydantic import BaseModel
from typing import List, Set
from rule import Rule 

class RuleData(BaseModel):
    num_undetermined_mines: int
    undetermined_cells: Set[str]
    
class SolverRequest(BaseModel):
    rules: List[RuleData]
    undetermined_mine_count: int
    num_informed_cells: int
    

app = FastAPI()

@app.post("/")
def get_frequenies(request: SolverRequest):
    
    # Helper to convert from the pydantic RuleData to the Rule class
    rules = [Rule(rule.num_undetermined_mines, frozenset(rule.undetermined_cells)) for rule in request.rules]
    
    # Determines the frequencies of the mines in the rules sent by the client
    ms = MinesweeperSolver(rules, request.undetermined_mine_count, request.num_informed_cells)
    frequencies, time = ms.solve()
     
    print(f"total time to solve: {time}")
    return frequencies