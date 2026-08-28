

export const saveToLs = (key, value) => {
    if (typeof window === 'undefined' || !window.localStorage) {return;}
    let storage = {};
    try {
        storage = JSON.parse(window.localStorage.getItem('dash_draggable')) || {};
    } catch (_error) {
        storage = {};
    }
    window.localStorage.setItem(
        `dash_draggable`,
        JSON.stringify({
            ...storage,
            [key]: value
        })
    );
};

export const getFromLs = (key) => {
    let ls = {};
    if (typeof window === 'undefined' || !window.localStorage) {return {};}
    try {
        ls = JSON.parse(window.localStorage.getItem(`dash_draggable`)) || {};
    } catch (_error) {
        /* Ignore */
    }
    return ls[key];
};
