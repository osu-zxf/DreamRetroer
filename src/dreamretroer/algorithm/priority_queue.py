import heapq

class PriorityQueue:
    def __init__(self):
        self._queue = []
        self._index = 0
        self._seen = {}

    def push(self, item, priority):
        if item in self._seen:
            self.remove(item)
        entry = [priority, self._index, item]
        self._index += 1
        self._seen[item] = entry
        heapq.heappush(self._queue, entry)

    def pop(self):
        if self.is_empty():
            return None
        else:
            [priority, index, item] = heapq.heappop(self._queue)
            while item == '<removed>':
                [priority, index, item] = heapq.heappop(self._queue)
            del self._seen[item]
            return priority, item

    def size(self):
        return len(self._seen)

    def is_empty(self):
        return len(self._seen) == 0

    def contains(self, item):
        return item in self._seen

    def remove(self, item):
        entry = self._seen.pop(item)
        entry[-1] = '<removed>'

    def update_priority(self, item, new_priority):
        self.remove(item)
        self.push(item, new_priority)

    def print(self):
        print('queue:', self._queue)
        print('dict:', self._seen)

if __name__ == "__main__":
    q = PriorityQueue()
    q.push('a', 1)
    q.push('b', 2)
    q.push('c', 3)
    q.push('d', 2)

    print(q.pop())  # a
    q.print()
    q.update_priority('d', 4)

    print(q.pop())  # d
    print(q.pop())  # b
    q.print()
    print(q.pop())  # c
    print(q.pop())  # None