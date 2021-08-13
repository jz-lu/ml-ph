class Base:
    def __init__(self, v1, v2):
        self.v1 = v1; self.v2 = v2
    def base_hello(self):
        print(f"Hello from Base!")

class Derived(Base):
    def __init__(self, v1, v2, v3):
        print("Hello from Derived!")
        print(f"Derived v1={v1}, v2={v2}, v3={v3}")
        super().__init__(v1, v2)
        self.base_hello()
    
Derived(10,20,30)